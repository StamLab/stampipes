import json
import os
import sys
import argparse
import logging
import requests
try:
    from concurrent.futures import ThreadPoolExecutor
except ImportError:
    from futures import ThreadPoolExecutor

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

STAMPIPES = os.getenv('STAMPIPES', '~/stampipes')

script_options = {
    "quiet": False,
    "debug": False,
    "base_api_url": None,
    "token": None,
    "outfile": os.path.join(os.getcwd(), "run.bash"),
    "sample_script_basename": "run.bash",
    "qsub_prefix": ".proc",
    "queue": "queue2",
    "dry_run": False,
    "no_mask": False,
    "bases_mask": None,
    "lane_ids": [],
    "flowcell_label": None,
    "tag_slug": None,
}

def parser_setup():

    parser = argparse.ArgumentParser()

    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true",
        help="Don't print info messages to standard out.")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true",
        help="Print all debug messages to standard out.")

    parser.add_argument("-a", "--api", dest="base_api_url",
        help="The base API url, if not the default live LIMS.")
    parser.add_argument("-t", "--token", dest="token",
        help="Your authentication token.  Required.")

    parser.add_argument("--script_template", dest="script_template",
        help="The script template to use.")

    parser.add_argument("-o", "--outfile", dest="outfile",
        help="Append commands to run this alignment to this file.")
    parser.add_argument("-b", "--sample-script-basename", dest="sample_script_basename",
        help="Name of the script that goes after the sample name.")
    parser.add_argument("--lane", dest="lane_ids", type=int, action="append",
        help="Lane ID")

    parser.add_argument("--flowcell_label", dest="flowcell_label", help="Flowcell Label")
    parser.add_argument("--tag", dest="tag", help="Lanes tagged by")

    parser.add_argument("--qsub-prefix", dest="qsub_prefix",
        help="Name of the qsub prefix in the qsub job name.  Use a . in front to make it non-cluttery.")
    parser.add_argument("--queue", dest="queue",
        help="SLURM queue for jobs.")

    parser.add_argument("-n", "--dry-run", dest="dry_run", action="store_true",
        help="Take no action, only print messages.")
    parser.add_argument("--no-mask", dest="no_mask", action="store_true",
        help="Don't use any barcode mask.")
    parser.add_argument("--bases_mask", dest="bases_mask",
        help="Set a bases mask.")

    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser


class ProcessSetUp(object):

    def __init__(self, args, api_url, token): 

        self.token = token
        self.api_url = api_url
        self.qsub_scriptname = args.sample_script_basename
        self.qsub_prefix = args.qsub_prefix
        self.queue = args.queue
        self.outfile = args.outfile
        self.script_template = args.script_template
        self.dry_run = args.dry_run
        self.no_mask = args.no_mask
        self.session = requests.Session()
        self.session.headers.update({'Authorization': "Token %s" % self.token})

        self.pool = ThreadPoolExecutor(max_workers=10)

    def api_single_result(self, url_addition=None, url=None):

        if url_addition:
           url = "%s/%s" % (self.api_url, url_addition)

        request = self.session.get(url)

        if request.ok:
            logging.debug(request.json())
            return request.json()
        else:
            logging.error("Could not get data from %s" % url)
            logging.error(request)
            return None

    def api_list_result(self, url_addition=None, url=None):

        more = True
        results = []

        if url_addition:
            url = "%s/%s" % (self.api_url, url_addition)

        while more:

            logging.debug("Fetching more results for query %s" % url)

            request = self.session.get(url)

            if not request.ok:
                logging.error(request)
                return None
            more_results = request.json()
            results.extend(more_results["results"])
            if more_results["next"]:
                url = more_results["next"]
            else:
                more = False

        return results

    def get_lane_process_info(self, lane_id):

        info = self.session.get("%s/flowcell_lane/%d/processing_information" % (self.api_url, lane_id))

        if info.ok:
            logging.debug(info.json())
            return info.json()
        else:
            logging.error("Could not find processing info for lane %d\n" % lane_id)
            logging.error(info)
            sys.exit(1)

    def get_process_template(self, process_template_id):

        if not process_template_id:
            logging.critical("No process template for alignment %d\n" % self.alignment_id)
            sys.exit(1)

        info = self.session.get("%s/process_template/%d" % (self.api_url, process_template_id))

        if info.ok:
            logging.debug(info.json())
            return info.json()
        else:
            logging.error("Could not find processing template for ID %d\n" % process_template_id)
            sys.exit(1)

    def setup_flowcell(self, flowcell_label):

        lanes = self.api_list_result("flowcell_lane?flowcell__label=%s" % flowcell_label)

        self.setup_lanes([lane["id"] for lane in lanes])

    def setup_tag(self, tag_slug):

        lane_tags = self.api_list_result("tagged_object?content_type=40&tag__slug=%s" % tag_slug)

        self.setup_lanes([lane_tag["object_id"] for lane_tag in lane_tags])

    def setup_lanes(self, lane_ids):
        self.pool.map(self.setup_lane, lane_ids)

    def setup_lane(self, lane_id):

        processing_info = self.get_lane_process_info(lane_id)

        self.create_script(processing_info)

    def add_script(self, script_file, lane_id, flowcell_label, sample_name):

        if not self.outfile:
            logging.debug("Writing script to stdout")
            outfile = sys.stdout
        else:
            logging.debug("Logging script to %s" % self.outfile)
            outfile = open(self.outfile, 'a')

        outfile.write("cd %s && " % os.path.dirname(script_file))
        outfile.write("qsub -N \"%s%s-%s-LANE#%d\" -q %s -cwd -V -S /bin/bash %s\n\n" % (self.qsub_prefix, sample_name, flowcell_label, lane_id, self.queue, script_file))

        outfile.close()

    def get_script_template(self):

        return open(self.script_template, 'r').read()

    def create_script(self, processing_info):

        lane = processing_info["libraries"][0]

        if not "directory" in lane:
            logging.critical("No directory for lane %d" % lane["id"])
            return False

        fastq_directory = lane["directory"]

        barcode = "NoIndex" if lane['barcode_index'] is None else lane['barcode_index']
        spreadsheet_name = "%s_%s_L00%d" % (lane['samplesheet_name'], barcode, lane['lane'])

        if not os.path.exists(fastq_directory):
            logging.critical("fastq directory %s does not exist, cannot continue" % fastq_directory)
            return False 

        script_file = os.path.join( fastq_directory, "%s-%s" % (spreadsheet_name, self.qsub_scriptname) )

        if self.dry_run:
            logging.info("Dry run, would have created: %s" % script_file)
            return True

        try:
            outfile = open(script_file, 'w')
        except FileNotFoundError:
            logging.critical("Could not create script file %s" % script_file)
            return False

        self.add_script(script_file, lane["id"], processing_info["flowcell"]["label"], spreadsheet_name)

        outfile.write("set -e -o pipefail\n")
        outfile.write("export SAMPLE_NAME=%s\n" % spreadsheet_name)
        outfile.write("export ASSAY=%s\n" % lane['assay'])
        outfile.write("export READLENGTH=%s\n" % processing_info['flowcell']['read_length'])
        if processing_info['flowcell']['paired_end']:
            outfile.write("export PAIRED=True\n")
        else:
            outfile.write("unset PAIRED\n")

        # Process with UMI if the barcode has one and this is a dual index
        # flowcell
        if lane['barcode1'] and lane['barcode1']['umi'] and processing_info['flowcell']['dual_index']:
            outfile.write("export UMI=True\n")
        else:
            outfile.write("unset UMI\n")

        outfile.write("export FLOWCELL_LANE_ID=%s\n" % lane['id'])
        outfile.write("export FASTQ_DIR=%s\n" % fastq_directory)
        outfile.write("export FLOWCELL=%s\n" % processing_info['flowcell']['label'])

        outfile.write("\n")
        outfile.write(self.get_script_template())
        outfile.close()


def main(args = sys.argv):
    """This is the main body of the program that by default uses the arguments
from the command line."""

    parser = parser_setup()
    poptions = parser.parse_args()

    if poptions.quiet:
        logging.basicConfig(level=logging.WARNING, format=log_format)
        logging.getLogger("requests").setLevel(logging.WARNING)
    elif poptions.debug:
        logging.basicConfig(level=logging.DEBUG, format=log_format)
    else:
        # Set up the logging levels
        logging.basicConfig(level=logging.INFO, format=log_format)
        logging.getLogger("requests").setLevel(logging.WARNING)

    if not poptions.base_api_url and "LIMS_API_URL" in os.environ:
        api_url = os.environ["LIMS_API_URL"]
    elif poptions.base_api_url:
        api_url = poptions.base_api_url
    else:
        logging.error("Could not find LIMS API URL.\n")
        sys.exit(1)

    if not poptions.token and "LIMS_API_TOKEN" in os.environ:
        token = os.environ["LIMS_API_TOKEN"]
    elif poptions.token:
        token = poptions.token
    else:
        logging.error("Could not find LIMS API TOKEN.\n")
        sys.exit(1)

    process = ProcessSetUp(poptions, api_url, token)

    process.setup_lanes(poptions.lane_ids)

    if poptions.flowcell_label:
        process.setup_flowcell(poptions.flowcell_label)

    if poptions.tag:
        process.setup_tag(poptions.tag)

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
