import json
import os
import sys
import argparse
import logging
import requests
import textwrap
from collections import OrderedDict
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
    "align_ids": [],
    "flowcell": None,
    "tag": None,
    "outfile": os.path.join(os.getcwd(), "run.bash"),
    "sample_script_basename": "run.bash",
    "qsub_queue": "queue0",
    "qsub_prefix": ".proc",
    "dry_run": False,
    "no_mask": False,
    "bases_mask": None,
    "redo_completed": False,
    "qsub_priority": -50,
    "auto_aggregate": False,
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

    parser.add_argument("--alignment", dest="align_ids", type=int, action="append",
        help="Run for this particular alignment.")
    parser.add_argument("--flowcell", dest="flowcell_label",
        help="Run for this particular flowcell label.")
    parser.add_argument("--tag", dest="tag",
        help="Run for alignments tagged here.")
    parser.add_argument("--project", dest="project",
        help="Run for alignments in this project.")

    parser.add_argument("--script_template", dest="script_template",
        help="The script template to use.")
    parser.add_argument("--qsub_priority", dest="qsub_priority", type=int,
        help="The priority to give scripts we are submitting.")

    parser.add_argument("-o", "--outfile", dest="outfile",
        help="Append commands to run this alignment to this file.")
    parser.add_argument("-b", "--sample-script-basename", dest="sample_script_basename",
        help="Name of the script that goes after the sample name.")
    parser.add_argument("--qsub-prefix", dest="qsub_prefix",
        help="Name of the qsub prefix in the qsub job name.  Use a . in front to make it non-cluttery.")
    parser.add_argument("--qsub-queue", dest="qsub_queue",
        help="Name of the SLURM partition")
    parser.add_argument("-n", "--dry-run", dest="dry_run", action="store_true",
        help="Take no action, only print messages.")
    parser.add_argument("--no-mask", dest="no_mask", action="store_true",
        help="Don't use any barcode mask.")
    parser.add_argument("--bases_mask", dest="bases_mask",
        help="Set a bases mask.")
    parser.add_argument("--redo_completed", dest="redo_completed", help="Redo alignments marked as completed.",
        action="store_true")
    parser.add_argument("--auto_aggregate", dest="auto_aggregate", help="Script created will also run auto-aggregations after alignments finished.",
        action="store_true")
    parser.add_argument("--align_base_dir", dest="align_base_dir", help="Create the alignment directory in this directory")

    parser.add_argument("--listout", dest="simple_output", help="Write only a list of alignments to run, rather than a script to submit them", action="store_true")

    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser


class ProcessSetUp(object):

    def __init__(self, args, api_url, token):

        self.token = token
        self.api_url = api_url
        self.qsub_scriptname = args.sample_script_basename
        self.qsub_prefix = args.qsub_prefix
        self.outfile = args.outfile
        self.dry_run = args.dry_run
        self.no_mask = args.no_mask
        self.redo_completed = args.redo_completed
        self.script_template = args.script_template
        self.qsub_priority = args.qsub_priority
        self.qsub_queue = args.qsub_queue
        self.auto_aggregate = args.auto_aggregate
        self.align_base_dir = args.align_base_dir

        self.simple_output = args.simple_output

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

    def get_align_process_info(self, alignment_id):

        process_info = self.api_single_result("flowcell_lane_alignment/%d/processing_information/" % alignment_id)

        if not process_info:
            logging.critical("Could not find processing info for alignment %d\n" % alignment_id)
            logging.critical(process_info)
            sys.exit(1)

        return process_info

    def get_process_template(self, align_id, process_template_id):

        if not process_template_id:
            logging.critical("No process template for alignment %d\n" % align_id)
            return None

        info = self.api_single_result("process_template/%d/" % (process_template_id))

        if not info:
            logging.critical("Could not find processing template for ID %d\n" % process_template_id)
            sys.exit(1)

        return info

    # Run alignment setup in parallel
    def setup_alignments(self, align_ids):
        for id, error in self.pool.map(self.setup_alignment, align_ids):
            if error:
                logging.debug("ALN%d result received, error: %s" % (id, error))
            else:
                logging.debug("ALN%d result received, OK" % id)

    def setup_alignment(self, align_id):

        try:
            processing_info = self.get_align_process_info(align_id)
            alignment = self.api_single_result("flowcell_lane_alignment/%d/" % (align_id))

            if self.redo_completed or not alignment['complete_time']:
                self.create_script(processing_info, alignment["id"])
                return (align_id, None)
            else:
                logging.info("Skipping completed alignment %d" % align_id)
                return (align_id, None)
        except Exception as e:
            logging.exception("Could not set up alignment %d}: (%s)" % (align_id, e))
            return (align_id, e)

    def get_lane_file(self, lane_id, purpose):
        candidates = self.api_list_result("file/?content_type=40&purpose__slug=%s&object_id=%d" % (purpose, lane_id))

        if not candidates:
            return None
        if len(candidates) > 1:
            return None

        return candidates[0]

    def setup_tag(self, tag_slug):

        align_tags = self.api_list_result("tagged_object/?content_type=47&tag__slug=%s" % tag_slug)

        self.setup_alignments([align_tag["object_id"] for align_tag in align_tags])

    def setup_project(self, project_id):
        logging.info("Setting up project #%s" % project_id)
        alignments = self.api_list_result("flowcell_lane_alignment/?lane__sample__project=%s" % project_id)
        self.setup_alignments([alignment["id"] for alignment in alignments])

    def setup_flowcell(self, flowcell_label):
        logging.info("Setting up flowcell for %s" % flowcell_label)
        alignments = self.api_list_result("flowcell_lane_alignment/?lane__flowcell__label=%s&page_size=1000" % flowcell_label)
        if self.auto_aggregate:
            for alignment in alignments:
                self.setup_alignment(alignment["id"])
            self.auto_aggregation_script(flowcell_label,alignments)
        else:
            self.setup_alignments([alignment["id"] for alignment in alignments])

    def auto_aggregation_script(self,flowcell_label,alignments):
        aaname_sentinel = "auto_agg_sentinel.%s" % (flowcell_label)

        if not self.outfile:
            logging.debug("Writing script to stdout")
            outfile = sys.stdout
        else:
            logging.debug("Logging script to %s" % self.outfile)
            outfile = open(self.outfile, 'a')

        contents = textwrap.dedent("""\
                   cd "$FLOWCELLS"/FC{label}_*
                   sentinel_dependencies=$(echo $PROCESSING | sed -e 's/,/,afterany:/g' | sed -e 's/^,afterany/--dependency=afterany/g')
                   sbatch --export=ALL -J {job_name} -o {job_name}.o%A -e {job_name}.e%A --partition={queue} --cpus-per-task=1 --ntasks=1 $sentinel_dependencies --mem-per-cpu=1000 --parsable --oversubscribe <<__AUTOAGG1__
                   #!/bin/bash
                   rm -f run_aggregations.bash
                   python $STAMPIPES/scripts/aggregateprocess.py --flowcell {label} --outfile run_aggregations.bash --qsub-queue {qqueue}
                   bash run_aggregations.bash
                   __AUTOAGG1__
                   """.format(label=flowcell_label,
                            job_name=aaname_sentinel,
                            queue=self.qsub_queue,
                            qqueue=self.qsub_queue))

        outfile.write(contents)
        outfile.close()

    def add_script(self, align_id, processing_info, script_file, sample_name):

        ram_megabytes = 2000

        if not self.outfile:
            logging.debug("Writing script to stdout")
            outfile = sys.stdout
        else:
            logging.debug("Logging script to %s" % self.outfile)
            outfile = open(self.outfile, 'a')

        if self.simple_output:
            outfile.write(script_file + "\n")
        else:
            outfile.write("cd %s && " % os.path.dirname(script_file))
            fullname = "%s%s-%s-ALIGN#%d" % (self.qsub_prefix,sample_name,processing_info['flowcell']['label'],align_id)
            outfile.write("jobid=$(sbatch --export=ALL -J %s -o %s.o%%A -e %s.e%%A --partition=%s --cpus-per-task=1 --ntasks=1 --mem-per-cpu=%d --parsable --oversubscribe <<__ALIGNPROC__\n#!/bin/bash\nbash %s\n__ALIGNPROC__\n)\nPROCESSING=\"$PROCESSING,$jobid\"\n\n" % (fullname, fullname, fullname, self.qsub_queue, ram_megabytes, script_file))
        outfile.close()

    def get_script_template(self, process_template):

        if self.script_template:
            script_path = self.script_template
        else:
            script_path = os.path.expandvars(process_template["process_version"]["script_location"])
        return open(script_path, 'r').read()

    def create_script(self, processing_info, align_id):

        lane = processing_info["libraries"][0]
        alignment = [a for a in lane["alignments"] if a["id"] == align_id][0]

        if not "process_template" in alignment:
            logging.error("Alignment %d has no process template" % align_id)
            return False

        process_template = self.get_process_template(align_id, alignment["process_template"])

        if not process_template:
            return False

        flowcell_directory = processing_info['flowcell']['directory']

        share_dir = lane.get("project_share_directory")
        if share_dir:
            flowcell_directory = os.path.join(share_dir, "alignments")
        if not flowcell_directory:
            logging.error("Alignment %d has no flowcell directory for flowcell %s" % (align_id, processing_info['flowcell']['label']))
            return False

        fastq_directory = os.path.join(flowcell_directory, "Project_%s" % lane['project'], "Sample_%s" % lane['samplesheet_name'])

        # Reset the alignment's sample name if we decied not to use the barcode index mask
        if self.no_mask:
            alignment['sample_name'] = "%s_%s_L00%d" % (lane['samplesheet_name'], lane['barcode_index'], lane['lane'])

        align_dir = "align_%d_%s_%s" % (alignment['id'], alignment['genome_index'], alignment['aligner'])
        if alignment['aligner_version']:
            align_dir = "%s-%s" % (align_dir, alignment['aligner_version'])

        script_directory = os.path.join(fastq_directory, align_dir)
        if self.align_base_dir:
            script_directory = os.path.join(self.align_base_dir, align_dir)

        r1_fastq = self.get_lane_file(lane["id"], "r1-fastq")

        if not r1_fastq:
            logging.error("Missing r1-fastq for lane %d (alignment %d) - check dir %s" % (lane["id"], alignment["id"], fastq_directory))
            return False

        if processing_info['flowcell']['paired_end']:
            r2_fastq = self.get_lane_file(lane["id"], "r2-fastq")
            if not r2_fastq:
                logging.error("Missing r2-fastq for lane %d (alignment %d)" % (lane["id"], alignment["id"]))
                return False

        script_file = os.path.join( script_directory, "%s-%s" % (alignment['sample_name'], self.qsub_scriptname) )
        logging.info("Will write to %s" % script_file)


        # Set up & add environment variables
        env_vars = OrderedDict()

        env_vars["SAMPLE_NAME"] = alignment['sample_name']
        env_vars["BWAINDEX"]    = alignment['genome_index_location']
        env_vars["GENOME"]      = alignment['genome_index']
        env_vars["ASSAY"]       = lane['assay']
        env_vars["READLENGTH"]  = processing_info['flowcell']['read_length']
        if processing_info['libraries'] and processing_info['libraries'][0] and processing_info['libraries'][0].get('library_kit_method'):
            env_vars["LIBRARY_KIT"] = '"' + processing_info['libraries'][0]['library_kit_method'] + '"'
        else:
            env_vars["LIBRARY_KIT"] = None

        if processing_info['flowcell']['paired_end']:
            env_vars["PAIRED"] = "True"
        else:
            env_vars["PAIRED"] = None

        env_vars["FLOWCELL_LANE_ID"] = lane['id']
        env_vars["ALIGNMENT_ID"]     = alignment['id']
        env_vars["ALIGN_DIR"]        = os.path.join(fastq_directory, align_dir)
        env_vars["R1_FASTQ"]         = r1_fastq["path"]

        if processing_info['flowcell']['paired_end']:
            env_vars["R2_FASTQ"] = r2_fastq["path"]

        env_vars["FASTQ_DIR"] = fastq_directory
        env_vars["FLOWCELL"]  = processing_info['flowcell']['label']

        if "barcode1" in lane and lane["barcode1"]:
            p7_adapter = lane['barcode1']['adapter7']
            p5_adapter = lane['barcode1']['adapter5']
            if "barcode2" in lane and lane['barcode2']:
                # Override the "default" end adapter from barcode1
                p5_adapter = lane['barcode2']['adapter5_reverse_complement']

            if not p7_adapter or not p5_adapter:
                logging.warn("Alignment %d missing adapters, some processes might not work" % alignment['id'])

            env_vars["ADAPTER_P7"] = p7_adapter
            env_vars["ADAPTER_P5"] = p5_adapter

            # Process with UMI if the barcode has one and this is a dual index
            # flowcell
            if lane['barcode1']['umi'] and processing_info['flowcell']['dual_index']:
                env_vars["UMI"] = "True"
            else:
                env_vars["UMI"] = None
            env_vars["UMI_METHOD"] = lane['barcode1']['umi_method']

        # Set process template env var overrides
        if 'process_variables' in process_template and process_template['process_variables']:
            try:
                process_template_variables = json.loads(process_template['process_variables'],
                                                        object_pairs_hook=OrderedDict)
                for var, value in process_template_variables.items():
                    env_vars[var] = value
            except ValueError as e:
                logging.error("Could not parse process variables for align %d (template %d): '%s'" %
                              (
                                  alignment['id'],
                                  process_template['id'],
                                  process_template['process_variables']
                              ))
                return False

        if self.dry_run:
            logging.info("Dry run, would have created: %s" % script_file)
            logging.debug(env_vars)
            return True

        if not os.path.exists(script_directory):
            logging.info("Creating directory %s" % script_directory)
            os.makedirs(script_directory)

        # Append to master script
        self.add_script(align_id, processing_info, script_file, alignment['sample_name'])

        # Write file
        outfile = open(script_file, 'w')
        outfile.write("set -e -o pipefail\n")

        # Set env vars
        for var, value in env_vars.items():
            if value is not None:
                outfile.write("export %s=%s\n" % (var, value))
            else:
                outfile.write("unset %s\n" % var)

        outfile.write("\n")
        outfile.write("export QUEUE=%s\n" % (self.qsub_queue))
        outfile.write("\n")
        outfile.write(self.get_script_template(process_template))
        outfile.close()


def main(args = sys.argv):
    """This is the main body of the program that by default uses the arguments
from the command line."""

    parser = parser_setup()
    poptions = parser.parse_args()

    if poptions.quiet:
        logging.basicConfig(level=logging.WARNING, format=log_format)
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

    process.setup_alignments(poptions.align_ids)

    if poptions.flowcell_label:
        process.setup_flowcell(poptions.flowcell_label)

    if poptions.tag:
        process.setup_tag(poptions.tag)

    if poptions.project:
        process.setup_project(poptions.project)


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
