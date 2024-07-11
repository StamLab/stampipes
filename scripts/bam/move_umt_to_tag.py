#!/usr/bin/env python3
"""
move_umt_to_tag.py
Takes a BAM file as input, and produces one as output.
Looks for a UMT tag embedded in the read name, and moves it to the 'RX' tag.
Preserves order, and is a no-op on reads without a UMT tag (no '#' in name)
"""

import argparse

import pysam


def move_umi(read, tag):
    """
    Looks for the UMI embeded in the read name, places it in a tag and trims
    the read name
    """
    umi_loc = read.query_name.find("#")
    if umi_loc > -1:
        read.set_tag(tag, read.query_name[umi_loc + 1 :])
        read.query_name = read.query_name[:umi_loc]
    return read


def main():
    """Sets up parsing and runs the program"""
    # Parsing
    parser = argparse.ArgumentParser(
        prog="move_umt_to_tag",
        description="Moves the UMT/UMI out of a read name and into a BAM tag",
    )
    parser.add_argument(
        "input_alignment", type=str, help="Input alignment file (with UMT in name)"
    )
    parser.add_argument(
        "output_alignment", type=str, help="Output alignment file (with UMT in tag)"
    )
    parser.add_argument(
        "--tagname", type=str, default="RX", help="Name of tag to store UMI info int"
    )
    args = parser.parse_args()

    input_alignment = pysam.AlignmentFile(args.input_alignment, "rb")
    output_alignment = pysam.AlignmentFile(
        args.output_alignment, "wb", template=input_alignment
    )
    reads = input_alignment.fetch(until_eof=True)

    # Do the work
    for read in reads:
        output_alignment.write(move_umi(read, args.tagname))

    # clean-up and close files
    input_alignment.close()
    output_alignment.close()


if __name__ == "__main__":
    main()
