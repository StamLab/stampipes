#!/usr/bin/env python3

import argparse
import datetime
import json
import re
import sys
import textwrap

from collections import defaultdict

# requires BioPython which seems to be in our environment
# but only to reverse complement which we could figure out
# another way to do

# Usage: $0 -p processing.json

SCRIPT_OPTIONS = {
    "processing": "processing.json",
    "reverse_barcode1": False,
    "reverse_barcode2": False,
    "filename": "SampleSheet.withmask.{mask}.csv",
}


def parser_setup():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--processing",
        dest="processing",
        help="The JSON file to read barcodes from (default: processing.json)",
    )
    parser.add_argument(
        "--reverse_barcode1",
        dest="reverse_barcode1",
        action="store_true",
        help="Use reverse sequence for barcode1",
    )
    parser.add_argument(
        "--reverse_barcode2",
        dest="reverse_barcode2",
        action="store_true",
        help="Use reverse sequence for barcode2",
    )
    parser.add_argument(
        "--filename",
        help="The template to use for filename, with the {mask} formatting",
    )
    parser.set_defaults(**SCRIPT_OPTIONS)
    return parser


def get_barcode_assignments(
    data: dict, reverse_barcode1: bool, reverse_barcode2: bool
) -> "[dict]":
    assignments = []

    for libdata in data["libraries"]:
        assignment = {
            "lane": libdata.get("lane"),
            "sample": libdata.get("samplesheet_name"),
            "barcode1": "",
            "barcode2": "",
        }
        if assignment["sample"] == "None":
            assignment["sample"] = "LANE%d" % libdata["id"]
        if libdata.get("barcode1") is not None:
            assignment["barcode1"] = (
                libdata["barcode1"]["reverse_sequence"]
                if reverse_barcode1
                else libdata["barcode1"]["sequence"]
            )
        if libdata.get("barcode2") is not None:
            assignment["barcode2"] = (
                libdata["barcode2"]["reverse_sequence"]
                if reverse_barcode2
                else libdata["barcode2"]["sequence"]
            )

        assignments.append(assignment)

    return assignments


def make_samplesheet_header(name: str, date: str) -> str:
    template = textwrap.dedent("""\
    [Header]
    Investigator Name,{name}
    Project Name,{name}
    Experiment Name,{name}
    Date,{date}
    Workflow,GenerateFASTQ

    [Settings]

    [Data]
    Lane,SampleID,SampleName,index,index2
    """)
    return template.format(date=date, name=name)


def group_assignments(assignments: "[dict]") -> "[[dict]]":
    """Groups the barcode assignments by length"""
    barcode_length_combinations = defaultdict(list)

    def get_len(d):
        return 0 if d is None else len(d)

    for assignment in assignments:
        key = (
            get_len(assignment["barcode1"]),
            get_len(assignment["barcode2"]),
        )
        barcode_length_combinations[key].append(assignment)
    return barcode_length_combinations


def parse_mask(mask: str) -> "[[(str, int)]]":
    parts = []
    str_parts = mask.split(",")
    regex = r"(?P<letter>[yni])(?P<num>[0-9]*)"
    for part in str_parts:
        pieces = []
        for match in re.finditer(regex, part, flags=re.I):
            letter = match.group("letter").lower()
            num_str = match.group("num")
            num = 1 if len(num_str) == 0 else int(num_str)
            if pieces and pieces[-1][0] == letter:
                # Collapse same-letter adjacent pieces
                pieces[-1][1] += num
            else:
                pieces.append((letter, num))
        parts.append(pieces)
    return parts


def mask_to_str(mask: "[[(str, int)]]") -> str:
    """Convert a mask in parts back into a string"""

    def format_piece(letter, num):
        if num == 0:
            return ""
        elif num == 1:
            return letter
        else:
            return letter + str(num)

    return ",".join(
        ["".join([format_piece(*piece) for piece in part]) for part in mask]
    )


def adjust_mask_for_lengths(mask_parts, len1, len2):
    """
    Takes in a barcode-mask (in parts) and the barcode length, and adjusts the
    values of 'i' to match.
    """
    new_mask = []
    index_reads_seen = 0
    for read in mask_parts:
        read_len = sum(piece[1] for piece in read)
        is_index_read = any(piece[0] == "i" for piece in read)
        if is_index_read:
            if any(piece[0] == "y" for piece in read):
                raise Exception(
                    "Mixed read/index in barcode mask '{}', don't know how to deal with this".format(
                        mask_to_str(mask_parts)
                    )
                )
            index_reads_seen += 1
            if index_reads_seen == 1:
                # first barcode
                if len1 == read_len:
                    new_mask.append([("i", len1)])
                elif len1 > read_len:
                    new_mask.append([("i", read_len)])
                elif len1 < read_len:
                    new_mask.append([("i", len1), ("n", read_len - len1)])
            elif index_reads_seen == 2:
                # second barcode
                if len2 == read_len:
                    new_mask.append([("i", len2)])
                elif len2 > read_len:
                    new_mask.append([("i", read_len)])
                elif len2 < read_len:
                    new_mask.append([("i", len2), ("n", read_len - len2)])
        else:
            new_mask.append(read)
    return new_mask


def write_samplesheets(name, filename_template, date, root_mask, assignments):
    """Write out the sample sheets"""
    mask_parts = parse_mask(root_mask)
    max_bclen1 = 0
    max_bclen2 = 0
    index_reads_seen = 0
    for read in mask_parts:
        read_len = sum(piece[1] for piece in read)
        is_index_read = any(piece[0] == "i" for piece in read)
        if is_index_read:
            index_reads_seen += 1
            if index_reads_seen == 1:
                max_bclen1 = read_len
            if index_reads_seen == 2:
                max_bclen2 = read_len

    for assign in assignments:
        assign["barcode1"] = assign["barcode1"][:max_bclen1]
        assign["barcode2"] = assign["barcode2"][:max_bclen2]
    # Trim barcodes to make sure they fit in the read

    groups = group_assignments(assignments)

    for barcode_lengths, assigns in groups.items():
        new_mask = adjust_mask_for_lengths(mask_parts, *barcode_lengths)
        header = make_samplesheet_header(name, date)
        body = make_samplesheet_body(assigns)
        samplesheet_contents = header + body
        filename = filename_template.format(mask=mask_to_str(new_mask))
        print(
            "Writing {filename} with {new_mask}".format(
                filename=filename, new_mask=mask_to_str(new_mask)
            )
        )
        with open(filename, "w") as f:
            f.write(samplesheet_contents)


def make_samplesheet_body(barcode_assignments: "[dict]") -> str:
    """Create samplesheet text from assignments"""
    lines = []
    for ba in barcode_assignments:
        line = ",".join(
            [
                str(ba["lane"]),
                ba["sample"],
                ba["sample"],
                str(ba["barcode1"]),
                str(ba["barcode2"]),
                "",
            ]
        )
        lines.append(line)
    return "\n".join(sorted(lines))


def main(args=sys.argv):
    parser = parser_setup()
    poptions = parser.parse_args()

    process_json = open(poptions.processing)
    data = json.load(process_json)
    process_json.close()

    assignments = get_barcode_assignments(
        data,
        poptions.reverse_barcode1,
        poptions.reverse_barcode2,
    )
    mask = data["alignment_group"]["bases_mask"]
    write_samplesheets(
        name="Altius",
        filename_template=poptions.filename,
        date=str(datetime.date.today()),
        root_mask=mask,
        assignments=assignments,
    )


if __name__ == "__main__":
    main()
