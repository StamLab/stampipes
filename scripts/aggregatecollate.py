import os
import sys
import argparse
import logging

sys.path.append("/home/audrakj/stamlims_api")
print(sys.path)
from stamlims_api import rest  # noqa: E402
from stamlims_api.lims import files  # noqa: E402


log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

STAMPIPES = os.getenv("STAMPIPES", "~/stampipes")

script_options = {
    "quiet": False,
    "debug": False,
    "base_api_url": None,
    "token": None,
    "aggregation_ids": [],
    "outfile": os.path.join(os.getcwd(), "aggregatecollate.bash"),
    "overwrite": False,
    "script_name": "fastq_collate.bash",
    "qsub_prefix": ".col",
    "dry_run": False,
    "aggregation_base_directory": None,
    "aggregation_directory": None,
    "script_template": os.path.join(
        STAMPIPES, "processes/fastq", "collate_aggregation_fastq.bash"
    ),
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
        "-a",
        "--api",
        dest="base_api_url",
        help="The base API url, if not the default live LIMS.  Required if not in the environment under LIMS_API_URL",
    )
    parser.add_argument(
        "-t",
        "--token",
        dest="token",
        help="Your authentication token.  Required if not in the environment under LIMS_API_TOKEN",
    )
    parser.add_argument(
        "--aggregation_base_directory",
        dest="aggregation_base_directory",
        help="The base directory to put aggregations in.  Can get from environment AGGREGATIONS variable",
    )

    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        help="Append commands to run this alignment to this file.",
    )
    parser.add_argument(
        "--overwrite",
        dest="overwrite",
        action="store_true",
        help="Create a new outfile instead of appending commands.",
    )
    parser.add_argument(
        "--script_template", dest="script_template", help="The script template to use."
    )
    parser.add_argument(
        "--aggregation_directory",
        dest="aggregation_directory",
        help="The directory for the aggregation.  Will deduce if not given.",
    )
    parser.add_argument(
        "-b",
        "--script_basename",
        dest="script_name",
        help="Name of the script that goes after the sample name.",
    )

    parser.add_argument("--tag", dest="tag", help="Run for alignments tagged here.")
    parser.add_argument(
        "--aggregation",
        dest="aggregation_ids",
        type=int,
        action="append",
        help="Run for these aggregations (can be used more than once).",
    )

    parser.add_argument(
        "--qsub-prefix",
        dest="qsub_prefix",
        help="Name of the qsub prefix in the qsub job name.  Use a . in front to make it non-cluttery.",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        dest="dry_run",
        action="store_true",
        help="Take no action, only print messages.",
    )

    parser.set_defaults(**script_options)
    parser.set_defaults(quiet=False, debug=False)

    return parser


class ProcessSetUp(object):
    def __init__(self, args, aggregation_base_directory):
        self.api = rest.setup_api()
        self.qsub_scriptname = args.script_name
        self.qsub_prefix = args.qsub_prefix
        self.outfile = args.outfile
        self.dry_run = args.dry_run
        self.aggregation_base_directory = aggregation_base_directory
        self.aggregation_directory = args.aggregation_directory
        self.script_template = args.script_template
        self.overwrite = args.overwrite

    def get_aggregation_info(self, aggregation_id):
        results = self.api.single_result("aggregation/%d" % aggregation_id)

        if not results:
            logging.error(
                "Could not find information for aggregation %d" % aggregation_id
            )
            return None

        return results

    def get_aggregation_lanes(self, aggregation_id):
        results = self.api.list_result(
            "aggregation_lane/?aggregation=%d&include=True" % aggregation_id
        )

        if not results:
            logging.error("Could not find lanes for aggregation %d" % aggregation_id)
            return []

        return results

    def get_lane_fastq_file(self, aggregation_id, lane_id, file_purpose):
        logging.info(
            "Fetching files for alignment %d (Aggregation %d)"
            % (lane_id, aggregation_id)
        )

        results = files.get_object_files(
            self.api,
            object_id=lane_id,
            object_content_type=40,
            purpose_shorthand=file_purpose,
        )

        if not results:
            logging.error("Improperly executed file query")
            return None

        if len(results) != 1:
            logging.error(
                "Found %d files for alignment %d, require 1 (Aggregation %d)"
                % (len(results), lane_id, aggregation_id)
            )
            logging.error(results)
            return None
        file_info = results[0]

        return (file_info["path"], file_info["md5sum"])

    def get_library_info(self, aggregation_info):
        library_info = self.api.single_result(url=aggregation_info["library"])
        if not library_info:
            logging.critical(
                "Cannot proceed without library!  Could not get info from %s (Aggregation %d)"
                % (aggregation_info["library"], aggregation_info["id"])
            )
            sys.exit(1)
        return library_info

    def get_script_template(self, script_template):
        logging.info("Using script template %s" % script_template)
        return open(script_template, "r").read()

    def get_example_flowcell(self, aggregation_id, aggregation_lanes):
        included = None
        for aggregation_lane in aggregation_lanes:
            if aggregation_lane["include"]:
                included = aggregation_lane
                break
        if included is None:
            logging.error(
                "No included aggregation lanes for aggregation %s", aggregation_id
            )
            lane = None
        else:
            lane = self.api.single_result(url=included["lane"])

        if not lane:
            logging.critical(
                "Was not able to fetch lane %s (Aggregation %d)"
                % (aggregation_lane["lane"], aggregation_id)
            )
            sys.exit(1)

        flowcell = self.api.single_result(url=lane["flowcell"])
        if not flowcell:
            logging.critical(
                "Could not get flowcell at %d (Aggregation %d)"
                % (lane["flowcell"], aggregation_id)
            )
            sys.exit(1)

        return flowcell

    def add_script(self, aggregation_id, aggregation_folder, library_number):
        if self.overwrite:
            mode = "w"
        else:
            mode = "a"

        with open(self.outfile, mode) as runfile:
            runfile.write("cd %s && " % aggregation_folder)
            runfile.write(
                'qsub -V -cwd -S /bin/bash -N "%sLN%d_AGG#%d" %s\n\n'
                % (
                    self.qsub_prefix,
                    library_number,
                    aggregation_id,
                    self.qsub_scriptname,
                )
            )

    def get_aggregation_directory(self, aggregation):
        return files.get_object_directories(
            self.api, object=aggregation, purpose_shorthand="aggregation-directory"
        )[0]["path"]

    def setup_tag(self, tag_slug):
        query_arguments = {"content_type": 126, "tag__slug": tag_slug}

        aggregation_tags = self.api.list_result(
            url_addition="tagged_object/", query_arguments=query_arguments
        )

        [
            self.setup_aggregation(aggregation_tag["object_id"])
            for aggregation_tag in aggregation_tags
        ]

    def setup_aggregation(self, aggregation_id):
        aggregation = self.get_aggregation_info(aggregation_id)

        if not aggregation:
            return False

        aggregation_lanes = self.get_aggregation_lanes(aggregation_id)

        if not aggregation_lanes:
            return False

        library_info = self.get_library_info(aggregation)
        aggregation_folder = self.get_aggregation_directory(aggregation)
        # flowcell = self.get_example_flowcell(aggregation_id, aggregation_lanes)

        logging.info("Aggregation %d folder: %s" % (aggregation_id, aggregation_folder))
        logging.debug(aggregation)

        missing = False

        r1_files = []
        r2_files = []
        for aggregation_lane in aggregation_lanes:
            lane_id = int(aggregation_lane["lane"].strip("/").split("/")[-1])
            if not aggregation_lane["include"]:
                logging.info(
                    "Not including lane %s (Aggregation %d)" % (lane_id, aggregation_id)
                )
                continue
            alignment_endpoint = aggregation_lane["alignment"]

            if not alignment_endpoint:
                logging.info(
                    "Not including lane %s because no alignment set (Aggregation %d)"
                    % (lane_id, aggregation_id)
                )

            # alignment_id = int(alignment_endpoint.strip("/").split("/")[-1])

            r1_fastq = self.get_lane_fastq_file(aggregation_id, lane_id, "r1-fastq")
            r2_fastq = self.get_lane_fastq_file(aggregation_id, lane_id, "r2-fastq")

            if not r1_fastq or not r2_fastq:
                logging.critical(
                    "Missing either R1: %s or R2: %s for alignment %s for lane %s, skipping (Aggregation %d)"
                    % (
                        str(r1_fastq),
                        str(r2_fastq),
                        alignment_endpoint,
                        lane_id,
                        aggregation_id,
                    )
                )
                missing = True
                continue
            else:
                logging.info(r1_fastq)
                logging.info(r2_fastq)
                r1_files.append(r1_fastq)
                r2_files.append(r2_fastq)

        if missing:
            return False

        script_contents = self.get_script_template(self.script_template)

        if not script_contents:
            logging.critical("No script contents")
            return

        os.makedirs(aggregation_folder, exist_ok=True)

        file_record = open("%s/r1.fastq.txt" % aggregation_folder, "w")
        file_record.write("\n".join(["\t".join(fastq) for fastq in r1_files]))
        file_record.close()

        file_record = open("%s/r2.fastq.txt" % aggregation_folder, "w")
        file_record.write("\n".join(["\t".join(fastq) for fastq in r2_files]))
        file_record.close()

        script_file = os.path.join(aggregation_folder, self.qsub_scriptname)
        logging.info("Creating script file %s" % script_file)

        script = open(script_file, "w")
        script.write("export AGGREGATION_ID=%d\n" % aggregation_id)
        script.write("export LIBRARY=%d\n" % library_info["number"])
        script.write("export LIBRARY_NAME=LN%d\n" % library_info["number"])
        script.write(
            'export R1_FILES="%s"\n' % " ".join([r1_file[0] for r1_file in r1_files])
        )
        script.write(
            'export R2_FILES="%s"\n' % " ".join([r2_file[0] for r2_file in r2_files])
        )
        script.write("export AGGREGATION_FOLDER=%s\n" % aggregation_folder)

        script.write(script_contents)

        script.close()

        self.add_script(aggregation_id, aggregation_folder, library_info["number"])


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

    if not poptions.aggregation_base_directory and "AGGREGATIONS" in os.environ:
        aggregation_base_dir = os.environ["AGGREGATIONS"]
    elif poptions.aggregation_base_directory:
        aggregation_base_dir = poptions.aggregation_base_directory
    else:
        aggregation_base_dir = None

    process = ProcessSetUp(poptions, aggregation_base_dir)

    for aggregation_id in poptions.aggregation_ids:
        process.setup_aggregation(aggregation_id)

    if poptions.tag:
        process.setup_tag(poptions.tag)


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
