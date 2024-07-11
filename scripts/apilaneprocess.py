#!/usr/bin/env python3
import argparse
import collections
import logging
import os
import re
import sys

try:
    from concurrent.futures import ThreadPoolExecutor
except ImportError:
    from futures import ThreadPoolExecutor

from stamlims_api import lims, rest

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

POOL_INFO = {}
SCRIPTS_WRITTEN = set()

STAMPIPES = os.getenv("STAMPIPES", "~/stampipes")

script_options = {
    "quiet": False,
    "debug": False,
    "outfile": os.path.join(os.getcwd(), "run.bash"),
    "sample_script_basename": "run.bash",
    "qsub_prefix": ".proc",
    "queue": "queue0",
    "dry_run": False,
    "no_mask": False,
    "bases_mask": None,
    "lane_ids": [],
    "flowcell_label": None,
    "tag_slug": None,
}


def parser_setup():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-q",
        "--quiet",
        dest="quiet",
        action="store_true",
        help="Don't print info messages to standard out.",
    )
    parser.add_argument(
        "-d",
        "--debug",
        dest="debug",
        action="store_true",
        help="Print all debug messages to standard out.",
    )

    parser.add_argument(
        "--script_template", dest="script_template", help="The script template to use."
    )

    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        help="Append commands to run this alignment to this file.",
    )
    parser.add_argument(
        "-b",
        "--sample-script-basename",
        dest="sample_script_basename",
        help="Name of the script that goes after the sample name.",
    )
    parser.add_argument(
        "--lane", dest="lane_ids", type=int, action="append", help="Lane ID"
    )

    parser.add_argument(
        "--flowcell_label", dest="flowcell_label", help="Flowcell Label"
    )
    parser.add_argument("--tag", dest="tag", help="Lanes tagged by")

    parser.add_argument(
        "--qsub-prefix",
        dest="qsub_prefix",
        help="Name of the qsub prefix in the qsub job name.  Use a . in front to make it non-cluttery.",
    )
    parser.add_argument("--queue", dest="queue", help="SLURM partition for jobs.")

    parser.add_argument(
        "-n",
        "--dry-run",
        dest="dry_run",
        action="store_true",
        help="Take no action, only print messages.",
    )
    parser.add_argument(
        "--no-mask",
        dest="no_mask",
        action="store_true",
        help="Don't use any barcode mask.",
    )
    parser.add_argument("--bases_mask", dest="bases_mask", help="Set a bases mask.")

    parser.set_defaults(**script_options)
    parser.set_defaults(quiet=False, debug=False)

    return parser


class ProcessSetUp(object):
    def __init__(self, args, api):
        self.api = api
        self.qsub_scriptname = args.sample_script_basename
        self.qsub_prefix = args.qsub_prefix
        self.queue = args.queue
        self.outfile = args.outfile
        self.script_template = args.script_template
        self.dry_run = args.dry_run
        self.no_mask = args.no_mask

        self.pool = ThreadPoolExecutor(max_workers=6)

    def get_lane_process_info(self, lane_id):
        info = self.api.get_single_result(
            url_addition="flowcell_lane/%d/processing_information" % (lane_id)
        )

        if info:
            return info
        else:
            logging.error("Could not find processing info for lane %d\n" % lane_id)
            sys.exit(1)

    def get_process_template(self, process_template_id):
        if not process_template_id:
            logging.critical(
                "No process template for alignment %d\n" % self.alignment_id
            )
            sys.exit(1)

        info = self.api.get_single_result(
            url_addition="process_template/%d" % (process_template_id)
        )

        if info:
            return info
        else:
            logging.error(
                "Could not find processing template for ID %d\n" % process_template_id
            )
            sys.exit(1)

    def setup_flowcell(self, flowcell_label):
        lanes = self.api.get_list_result(
            url_addition="flowcell_lane/",
            query_arguments={"flowcell__label": flowcell_label},
            page_size=1000,
            item_limit=50000,
        )

        if not lanes:
            logging.error("Flowcell %s has no lanes" % flowcell_label)
            return

        logging.debug(
            "Setting up flowcell %s with %d lanes" % (flowcell_label, len(lanes))
        )

        self.setup_lanes([lane["id"] for lane in lanes])

    def setup_tag(self, tag_slug):
        flowcelllane_contenttype = lims.content_types.contenttype_from_model_name(
            self.api, model_name="FlowcellLane"
        )
        lane_tags = self.api.get_list_result(
            url_addition="tagged_object",
            query_arguments={
                "content_type": flowcelllane_contenttype["id"],
                "tag__slug": tag_slug,
            },
        )

        if not lane_tags:
            logging.error("Tag %s has no lanes" % lane_tags)

        logging.debug("Setting up tag %s " % tag_slug)

        self.setup_lanes([lane_tag["object_id"] for lane_tag in lane_tags])

    def setup_lanes(self, lane_ids):
        logging.debug("Setting up lane IDs %s" % str(lane_ids))

        if len(lane_ids) != len(set(lane_ids)):
            logging.warning(
                "Duplicate lane IDs! %s "
                % [
                    item
                    for item, count in collections.Counter(lane_ids).items()
                    if count > 1
                ]
            )

        # self.pool.map(self.setup_lane, lane_ids)
        for lane_id in lane_ids:
            self.setup_lane(lane_id)

    def setup_lane(self, lane_id):
        logging.debug("Setting up lane %d" % lane_id)

        processing_info = self.get_lane_process_info(lane_id)

        pool_name = None

        if len(processing_info.get("libraries", [])) == 1 and processing_info[
            "libraries"
        ][0].get("samplesheet_name", "").startswith("LP"):
            pool_name = processing_info["libraries"][0]["samplesheet_name"]
            pool_number = int(pool_name[2:])  # remove leading LP
            pool_data = self.api.get_single_result(
                url_addition="library_pool/?number=%d" % pool_number
            )["results"][0]
            pool_id = pool_data["id"]
        else:
            try:
                lib_number = processing_info["libraries"][0]["library"]
                library_info = self.api.get_single_result(
                    url_addition="library/?number=%d" % lib_number
                )["results"][0]
                logging.debug("Info is %s", library_info)
                pools = library_info["librarypools"]
                if pools:
                    pool_name = pools[0]["object_name"]
                    pool_id = pools[0]["id"]
                    logging.debug("Lane %d is pool %s", lib_number, pool_name)
                else:
                    logging.debug("Lane %d is not pool", lib_number)
            except Exception:
                pass

        global POOL_INFO
        if pool_name and pool_name not in POOL_INFO:
            pool_data = self.api.get_single_result(
                url_addition="library_pool/%d/" % pool_id
            )
            bc1 = None
            bc2 = None
            if pool_data["barcode1"]:
                bc1 = self.api.get_single_result(url=pool_data["barcode1"])[
                    "reverse_sequence"
                ]
            if pool_data["barcode2"]:
                bc2 = self.api.get_single_result(url=pool_data["barcode2"])[
                    "reverse_sequence"
                ]
            barcode = "-".join(bc for bc in [bc1, bc2] if bc)
            POOL_INFO[pool_name] = {"barcode": barcode}

        self.create_script(processing_info, pool_name)

    def add_script(self, script_file, lane_id, flowcell_label, sample_name):
        # Hacks to deduplicate files written for library pools
        global SCRIPTS_WRITTEN
        if script_file in SCRIPTS_WRITTEN:
            return
        SCRIPTS_WRITTEN.add(script_file)

        if not self.outfile:
            logging.debug("Writing script to stdout")
            outfile = sys.stdout
        else:
            logging.debug("Logging script to %s" % self.outfile)
            outfile = open(self.outfile, "a")

        outfile.write("cd %s && " % os.path.dirname(script_file))
        fullname = "%s%s-%s-Lane#%d" % (
            self.qsub_prefix,
            sample_name,
            flowcell_label,
            lane_id,
        )
        outfile.write(
            "sbatch --export=ALL -J %s -o %s.o%%A -e %s.e%%A --partition=%s --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__LANEPROC__\n#!/bin/bash\nbash %s\n__LANEPROC__\n\n"
            % (fullname, fullname, fullname, self.queue, script_file)
        )
        outfile.close()

    def get_script_template(self):
        return open(self.script_template, "r").read()

    def create_script(self, processing_info, pool=None):
        lane = processing_info["libraries"][0]

        if "directory" not in lane:
            logging.critical("No directory for lane %d" % lane["id"])
            return False
        fastq_directory = lane["directory"]
        alt_dir = lane.get("project_share_directory", "")
        if alt_dir:
            fastq_directory = os.path.join(
                alt_dir,
                "fastq",
                "Project_%s" % lane["project"],
                "Sample_%s" % lane["samplesheet_name"],
            )

        if pool:
            flowcell_dir = re.sub(r"/Project.*", "", lane["directory"])
            if alt_dir:
                flowcell_dir = alt_dir
            fastq_directory = os.path.join(
                flowcell_dir, "Project_%s" % lane["project"], "LibraryPool_%s" % pool
            )

        barcode = "NoIndex" if lane["barcode_index"] is None else lane["barcode_index"]
        try:
            # Preferred name
            spreadsheet_name = lane["alignments"][0]["sample_name"]
        except (KeyError, IndexError):
            # Fallback method, doesn't always have the same barcode string
            spreadsheet_name = "%s_%s_L00%d" % (
                lane["samplesheet_name"],
                barcode,
                lane["lane"],
            )
            logging.warning(
                "No alignment sample_name for lane, using %s instead" % spreadsheet_name
            )

        if pool:
            global POOL_INFO
            barcode = POOL_INFO[pool]["barcode"]
            spreadsheet_name = "%s_%s_L00%d" % (pool, barcode, lane["lane"])
        # print("DBG:", pool, spreadsheet_name, POOL_INFO)

        if not os.path.exists(fastq_directory):
            logging.critical(
                "fastq directory %s does not exist, cannot continue" % fastq_directory
            )
            return False

        script_file = os.path.join(
            fastq_directory, "%s-%s" % (spreadsheet_name, self.qsub_scriptname)
        )

        if self.dry_run:
            logging.info("Dry run, would have created: %s" % script_file)
            return True

        try:
            outfile = open(script_file, "w")
        except FileNotFoundError:
            logging.critical("Could not create script file %s" % script_file)
            return False

        self.add_script(
            script_file,
            lane["id"],
            processing_info["flowcell"]["label"],
            spreadsheet_name,
        )

        outfile.write("set -e -o pipefail\n")
        outfile.write("export SAMPLE_NAME=%s\n" % spreadsheet_name)
        outfile.write("export ASSAY=%s\n" % lane["assay"])
        outfile.write(
            "export READLENGTH=%s\n" % processing_info["flowcell"]["read_length"]
        )
        if processing_info["flowcell"]["paired_end"]:
            outfile.write("export PAIRED=True\n")
        else:
            outfile.write("unset PAIRED\n")

        # Process with UMI if the barcode has one and this is a dual index
        # flowcell
        if (
            lane["barcode1"]
            and lane["barcode1"]["umi"]
            and processing_info["flowcell"]["dual_index"]
        ):
            outfile.write("export UMI=True\n")
        else:
            outfile.write("unset UMI\n")

        outfile.write("export FLOWCELL_LANE_ID=%s\n" % lane["id"])
        outfile.write("export FASTQ_DIR=%s\n" % fastq_directory)
        outfile.write("export FLOWCELL=%s\n" % processing_info["flowcell"]["label"])

        outfile.write("\n")
        outfile.write(self.get_script_template())
        outfile.close()


def main(args=sys.argv):
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

    api = rest.setup_api()

    process = ProcessSetUp(poptions, api)

    if len(poptions.lane_ids) > 0:
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
