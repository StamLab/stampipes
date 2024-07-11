#!/usr/bin/env python3
from __future__ import unicode_literals

import argparse
import glob
import json
import logging
import os
import re
from collections import defaultdict

LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

SCRIPT_OPTIONS = {
    "quiet": False,
    "debug": False,
    "base_dir": os.getcwd(),
    "processing_file": os.path.join(os.getcwd(), "processing.json"),
    "dry_run": False,
}


def parser_setup():
    """Sets up parser"""

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
        "-i", "--input-dir", dest="input_dir", help="The input directory to use."
    )
    parser.add_argument(
        "-o", "--output-dir", dest="output_dir", help="The output directory to use."
    )
    parser.add_argument(
        "-p",
        "--processing_file",
        dest="processing_file",
        help="The processing_file to use as a guide.",
    )
    parser.add_argument(
        "--merge-across-lanes",
        action="store_true",
        help="Merge across physical flowcell lanes. (Recommended for nextseq, not for Novaseq)",
    )

    parser.add_argument(
        "--dry-run",
        dest="dry_run",
        action="store_true",
        help="Only print out planned symlinks.",
    )

    parser.set_defaults(**SCRIPT_OPTIONS)
    parser.set_defaults(quiet=False, debug=False)

    return parser


def create_links(
    lane,
    read,
    input_basedir,
    output_basedir,
    dry_run=False,
    undetermined=False,
    is_pool=False,
    merge_across_lanes=False,
):
    """
    Create the links between the input directories and output dir
    If dry_run is passed, will print them instead of creating them
    """

    short_name = lane["samplesheet_name"]
    if lane.get("alignments"):
        sample_name = lane["alignments"][0]["sample_name"]
    else:
        bc1 = lane["barcode1"]["reverse_sequence"] if lane.get("barcode1") else ""
        bc2 = lane["barcode2"]["reverse_sequence"] if lane.get("barcode2") else ""
        lane_num = int(lane["lane"])
        if bc1 and bc2:
            sample_name = "%s_%s-%s_L%03d" % (short_name, bc1, bc2, lane_num)
        elif bc1:
            sample_name = "%s_%s_L%03d" % (short_name, bc1, lane_num)
        elif bc2:
            sample_name = "%s_%s_L%03d" % (short_name, bc2, lane_num)
        else:
            sample_name = "%s_L%03d" % (short_name, lane_num)

    if lane.get("library_pool"):
        is_pool = True

    if undetermined:
        output_dir = os.path.join(
            output_basedir, "Undetermined_indices", "Sample_lane1"
        )
    else:
        prefix = "LibraryPool" if is_pool else "Sample"
        output_dir = os.path.join(
            output_basedir,
            "Project_%s" % lane["project"],
            "%s_%s" % (prefix, lane["samplesheet_name"]),
        )

    short_name = re.sub(r"_", "-", short_name)
    lane_lane = (
        "*"
        if (merge_across_lanes or "lane" not in lane)
        else "_L%03d" % int(lane["lane"])
    )
    input_wildcard = os.path.join(
        input_basedir, "%s_S*%s_%s_???.fastq.gz" % (short_name, lane_lane, read)
    )

    if not dry_run and not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # This will fail if we have the same sample listed multiple times in the
    # samplesheet (run with different barcodes).
    # (Since we use the library letter in the sample name, this would imply
    # that the library has multiple barcodes, which isn't a thing that makes
    # sense in our system)
    input_fastq = sorted(glob.glob(input_wildcard))

    logging.debug("Looking for %s", input_wildcard)
    for idx, input_file in enumerate(input_fastq, start=1):
        output_name = "%s_%s_%03d.fastq.gz" % (sample_name, read, idx)
        output_file = os.path.join(output_dir, output_name)

        rel_path = os.path.relpath(input_file, output_dir)

        logging.info("Linking %s => %s" % (rel_path, output_file))
        if not dry_run and not os.path.exists(output_file):
            os.symlink(rel_path, output_file)


def main():
    """This is the main body of the program that by default uses the arguments
    from the command line."""

    parser = parser_setup()
    poptions = parser.parse_args()

    if poptions.quiet:
        logging.basicConfig(level=logging.WARNING, format=LOG_FORMAT)
    elif poptions.debug:
        logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)
    else:
        # Set up the logging levels
        logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)

    input_dir = poptions.input_dir

    data = json.loads(open(poptions.processing_file, "r").read())

    for lane in data["libraries"]:
        create_links(
            lane,
            "R1",
            input_dir,
            poptions.output_dir,
            poptions.dry_run,
            merge_across_lanes=poptions.merge_across_lanes,
        )
        create_links(
            lane,
            "R2",
            input_dir,
            poptions.output_dir,
            poptions.dry_run,
            merge_across_lanes=poptions.merge_across_lanes,
        )

    undet_lane = {
        "alignments": [{"sample_name": "lane1_Undetermined_L001"}],
        "samplesheet_name": "Undetermined",
    }
    for read in ["R1", "R2"]:
        create_links(
            undet_lane,
            read,
            input_dir,
            poptions.output_dir,
            poptions.dry_run,
            undetermined=True,
            merge_across_lanes=poptions.merge_across_lanes,
        )

    # Set up conversion table
    libs_to_lanes = defaultdict(set)
    for lane in data["libraries"]:
        libs_to_lanes[lane["library"]].add(lane["lane"])

    for pool, info in data["library_pools"].items():
        barcode = info["barcode1"]
        if info.get("barcode2"):
            barcode = "%s_%s" % (barcode, info["barcode2"])
        lane_nums = set()
        for lib in info["libraries"]:
            lib_num = int(re.sub(r"[^\d]+", "", lib))
            lane_nums.update(libs_to_lanes[lib_num])

        for lane_num in sorted(lane_nums):
            out_name = "%s_%s_L00%d" % (pool, barcode, lane_num)
            lane = {
                "samplesheet_name": pool,
                "alignments": [{"sample_name": out_name}],
                "project": "Lab",
            }
            for read in ["R1", "R2"]:
                create_links(
                    lane,
                    read,
                    input_dir,
                    poptions.output_dir,
                    poptions.dry_run,
                    is_pool=True,
                    merge_across_lanes=poptions.merge_across_lanes,
                )


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
