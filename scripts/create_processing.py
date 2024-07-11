#!/usr/bin/env python3
import argparse
import json
import logging
import os
import sys

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

STAMPIPES = os.getenv("STAMPIPES", "~/stampipes")

script_files = {
    "bwa": "%s/processes/bwa/process_bwa_paired_trimmed.bash" % STAMPIPES,
    "tophat": "%s/processes/tophat/process_tophat_paired.bash" % STAMPIPES,
    "fastqc": "%s/processes/fastqc_only.bash" % STAMPIPES,
}

# Right now, these are always run.
flowcell_script_files = {
    "barcodetally": "%s/processes/barcode_tally.bash" % STAMPIPES,
}

script_contents = {}

script_options = {
    "quiet": False,
    "debug": False,
    "process_config": "processing.json",
    "outfile": os.path.join(os.getcwd(), "run.bash"),
    "sample_script_basename": "run.bash",
    "qsub_prefix": ".proc",
    "template_script": None,
    "project_filter": [],
    "no_mask": False,
    "dry_run": False,
    "add_flowcell_scripts": False,
    "umi_filter": None,
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
        "-o",
        "--outfile",
        dest="outfile",
        help="The master script to run all sample scripts.",
    )
    parser.add_argument(
        "-p",
        "--process-config",
        dest="process_config",
        help="The process config to work off of.",
    )
    parser.add_argument(
        "-b",
        "--sample-script-basename",
        dest="sample_script_basename",
        help="Name of the script that goes after the sample name.",
    )
    parser.add_argument(
        "--project",
        dest="project_filter",
        action="append",
        help="Run for this particular project. Can be specified multiple times.",
    )

    parser.add_argument(
        "--sample",
        dest="sample_filter",
        type=int,
        action="append",
        help="Run for this particular sample. Can be specified multiple times.",
    )
    parser.add_argument(
        "--library",
        dest="library_filter",
        type=int,
        action="append",
        help="Run for this particular library. Can be specified multiple times.",
    )
    parser.add_argument(
        "--alignment",
        dest="alignment_filter",
        type=int,
        action="append",
        help="Run for this particular alignment. Can be specified multiple times.",
    )
    parser.add_argument(
        "--lane",
        dest="lanes",
        type=int,
        action="append",
        help="Run for these given lanes.  Can be specified multiple times.",
    )
    parser.add_argument(
        "--no_umi",
        dest="umi_filter",
        action="store_false",
        help="Run and don't include any UMI lanes.",
    )
    parser.add_argument(
        "--only_umi",
        dest="umi_filter",
        action="store_true",
        help="Run and only include UMI lanes.",
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
    parser.add_argument(
        "-t",
        "--template-script",
        dest="template_script",
        help="Template script to make for each valid library if not defaults",
    )
    parser.add_argument(
        "--no-mask",
        dest="no_mask",
        action="store_true",
        help="If this is set to true, remake SAMPLE_NAME with no barcode mask.",
    )
    parser.add_argument(
        "--add_flowcell_scripts",
        dest="add_flowcell_scripts",
        action="store_true",
        help="If this is set to true, add on flowcell level scripts to the end of the run collection.",
    )

    parser.add_argument(
        "-i",
        "--ignore_failed_lanes",
        dest="ignore_failed_lanes",
        action="store_true",
        help="Ignore failed lanes",
    )

    parser.set_defaults(**script_options)
    parser.set_defaults(quiet=False, debug=False)

    return parser


class ProcessSetUp(object):
    def __init__(self, args):
        self.processing_configfile = args.process_config
        self.qsub_scriptname = args.sample_script_basename
        self.qsub_prefix = args.qsub_prefix
        self.outfile = args.outfile
        self.template_script = args.template_script
        self.project_filter = args.project_filter
        self.library_filter = args.library_filter
        self.sample_filter = args.sample_filter
        self.alignment_filter = args.alignment_filter
        self.no_mask = args.no_mask
        self.dry_run = args.dry_run
        self.add_flowcell_scripts = args.add_flowcell_scripts
        self.ignore_failed_lanes = args.ignore_failed_lanes
        self.umi_filter = args.umi_filter
        self.filter_lanes = args.lanes

        if self.template_script:
            self.template_script_content = open(self.template_script, "r").read()

    def include_lane(self, lane):
        if self.umi_filter is not None:
            if lane["barcode1"] and lane["barcode1"]["umi"]:
                umi = True
            else:
                umi = False

            if self.umi_filter and not umi:
                return False
            if not self.umi_filter and umi:
                return False

        if self.filter_lanes and lane["lane"] not in self.filter_lanes:
            return False

        if self.ignore_failed_lanes and lane["failed"]:
            logging.debug(
                "Skipping %s, failed and we are ignoring failed lanes",
                lane["samplesheet_name"],
            )
            return False

        if self.project_filter and lane["project"] not in self.project_filter:
            logging.debug(
                "Skipping %s, not in project filter", lane["samplesheet_name"]
            )
            return False

        if self.library_filter and lane["library"] not in self.library_filter:
            logging.debug(
                "Skipping %s, not in library filter", lane["samplesheet_name"]
            )
            return False

        if self.sample_filter and lane["sample"] not in self.sample_filter:
            logging.debug("Skipping %s, not in sample filter", lane["samplesheet_name"])
            return False

        if (
            self.alignment_filter
            and lane["alignments"]
            and lane["alignments"][0]["id"] not in self.alignment_filter
        ):
            logging.debug(
                "Skipping %s, not in alignment filter", lane["samplesheet_name"]
            )
            return False

        return True

    def create(self):
        self.processing_scripts = dict()
        self.p = json.loads(open(self.processing_configfile, "r").read())

        for lane in self.p["libraries"]:
            if self.include_lane(lane):
                self.create_script(lane)

        if self.add_flowcell_scripts:
            logging.info("Adding flowcell level scripts.")
            for script in flowcell_script_files.values():
                self.create_flowcell_script(script)

        self.run_scripts()

    def add_script(self, script_file, sample_name, priority):
        if priority not in self.processing_scripts:
            self.processing_scripts[priority] = list()

        self.processing_scripts[priority].append((sample_name, script_file))

    def run_scripts(self):
        if self.dry_run:
            return

        outfile = open(self.outfile, "w")

        for priority in sorted(self.processing_scripts.keys(), reverse=True):
            outfile.write("# Priority %s\n" % str(priority))
            for sample_name, script_file in self.processing_scripts[priority]:
                outfile.write("cd %s && " % os.path.dirname(script_file))
                outfile.write(
                    "qsub -N %s%s-%s -cwd -V -S /bin/bash %s\n\n"
                    % (
                        self.qsub_prefix,
                        sample_name,
                        self.p["flowcell"]["label"],
                        script_file,
                    )
                )

        outfile.close()

    def get_script_template(self, lane):
        if self.template_script:
            return self.template_script_content

        alignment = lane["alignments"][0]
        if not alignment["aligner"]:
            logging.info("# FastQC only %s", lane["sample"])
            base_script = "fastqc"
        else:
            base_script = alignment["aligner"]
            logging.info("# Aligning %s with %s", lane["sample"], base_script)

        if base_script not in script_contents:
            script_contents[base_script] = open(script_files[base_script], "r").read()

        return script_contents[base_script]

    # Probably ripe for a refactoring
    def create_flowcell_script(self, inscript):
        script_directory = os.path.join(
            self.p["alignment_group"]["directory"], "flowcell_scripts"
        )

        if not os.path.exists(script_directory):
            logging.info("Creating directory %s", script_directory)
            os.makedirs(script_directory)

        script_file = os.path.join(script_directory, os.path.basename(inscript))

        # TODO: Figure out the appropriate priority here.
        self.add_script(script_file, "flowcell_script", 0)

        if self.dry_run:
            return

        outfile = open(script_file, "w")
        outfile.write("set -e -o pipefail\n")
        outfile.write("export READLENGTH=%s\n" % self.p["flowcell"]["read_length"])
        if self.p["flowcell"]["paired_end"]:
            outfile.write("export PAIRED=True\n")

        outfile.write("export FLOWCELL=%s\n" % self.p["flowcell"]["label"])
        outfile.write(
            "export FLOWCELL_DIR=%s\n" % self.p["alignment_group"]["directory"]
        )
        outfile.write(
            "export PROCESSING=%s\n" % os.path.abspath(self.processing_configfile)
        )
        outfile.write("\n")
        outfile.close()

        os.system("cat %s >> %s" % (inscript, script_file))

    def create_script(self, lane):
        if not lane["alignments"]:
            return False

        alignment = lane["alignments"][0]

        fastq_directory = os.path.join(
            self.p["alignment_group"]["directory"],
            "Project_%s" % lane["project"],
            "Sample_%s" % lane["samplesheet_name"],
        )

        # Reset the alignment's sample name if we decied not to use the barcode index mask
        if self.no_mask:
            alignment["sample_name"] = "%s_%s_L00%d" % (
                lane["samplesheet_name"],
                lane["barcode_index"],
                lane["lane"],
            )

        align_dir = "align_%d_%s_%s" % (
            alignment["id"],
            alignment["genome_index"],
            alignment["aligner"],
        )
        if alignment["aligner_version"]:
            align_dir = "%s-%s" % (align_dir, alignment["aligner_version"])

        script_directory = "%s/%s" % (fastq_directory, align_dir)

        if not os.path.exists(script_directory):
            logging.info("Creating directory %s", script_directory)
            os.makedirs(script_directory)

        script_file = os.path.join(
            script_directory, "%s-%s" % (alignment["sample_name"], self.qsub_scriptname)
        )
        logging.info(script_file)

        self.add_script(script_file, alignment["sample_name"], alignment["priority"])

        if self.dry_run:
            return

        outfile = open(script_file, "w")
        outfile.write("set -e -o pipefail\n")
        outfile.write("export SAMPLE_NAME=%s\n" % alignment["sample_name"])
        outfile.write("export BWAINDEX=%s\n" % alignment["genome_index_location"])
        outfile.write("export GENOME=%s\n" % alignment["genome_index"])
        outfile.write("export ASSAY=%s\n" % lane["assay"])
        outfile.write("export READLENGTH=%s\n" % self.p["flowcell"]["read_length"])
        if self.p["flowcell"]["paired_end"]:
            outfile.write("export PAIRED=True\n")
        outfile.write("export FLOWCELL_LANE_ID=%s\n" % lane["id"])
        outfile.write("export ALIGNMENT_ID=%s\n" % alignment["id"])
        outfile.write("export ALIGN_DIR=%s/%s\n" % (fastq_directory, align_dir))
        outfile.write("export FASTQ_DIR=%s\n" % fastq_directory)
        outfile.write("export FLOWCELL=%s\n" % self.p["flowcell"]["label"])
        if "barcode1" in lane and lane["barcode1"]:
            p7_adapter = lane["barcode1"]["adapter7"]
            p5_adapter = lane["barcode1"]["adapter5"]
            if "barcode2" in lane and lane["barcode2"]:
                # Override the "default" end adapter from barcode1
                # TODO: Make sure we want adapter7, double-check lims methods
                p5_adapter = lane["barcode2"]["adapter7"]
            outfile.write("export ADAPTER_P7=%s\n" % p7_adapter)
            outfile.write("export ADAPTER_P5=%s\n" % p5_adapter)

            # Process with UMI if the barcode has one and this is a dual index
            # flowcell
            if lane["barcode1"]["umi"] and self.p["flowcell"]["dual_index"]:
                outfile.write("export UMI=True\n")

        outfile.write("\n")
        outfile.write(self.get_script_template(lane))
        outfile.close()


def main(args=sys.argv):
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

    process = ProcessSetUp(poptions)
    process.create()


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
