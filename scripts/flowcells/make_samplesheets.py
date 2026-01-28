#!/usr/bin/env python3
"""
Generate bcl-convert compatible sample sheets from processing.json

This script creates sample sheets with the bcl-convert v2 format, including:
- OverrideCycles in [Settings] section (replaces bcl2fastq's --use-bases-mask)
- BarcodeMismatchesIndex1/2 in [Settings] section
- Proper [Data] section with Lane, Sample_ID, index, index2 columns

Usage: $0 -p processing.json [--mismatches 1,1] [--reverse_barcode1] [--reverse_barcode2] [--filename SampleSheet.withmask.{mask}.csv]
"""

import argparse
import datetime
import json
import re
import sys
import textwrap
from collections import defaultdict
from dataclasses import dataclass, field

SCRIPT_OPTIONS = {
    "processing": "processing.json",
    "reverse_barcode1": False,
    "reverse_barcode2": False,
    "filename": "SampleSheet.withmask.{mask}.csv",
    "mismatches": "1,1",
}


@dataclass
class BasesMask:
    """
    Represents a bases mask / cycle override for Illumina sequencing.

    A mask consists of multiple "reads" (segments), each containing pieces
    that describe how to handle bases:
    - 'y' = use for sequencing read
    - 'i' = use for index/barcode
    - 'n' = skip/ignore
    - 'u' = use for UMI

    Examples:
        "y151,i8,i8,y151" - paired-end 151bp with dual 8bp indexes
        "y76,i8,y76" - paired-end 76bp with single 8bp index
    """

    # Each read is a list of (letter, count) tuples
    # e.g., [[('y', 151)], [('i', 8)], [('i', 8)], [('y', 151)]]
    reads: list[list[tuple[str, int]]] = field(default_factory=list)

    @classmethod
    def parse(cls, mask_str: str) -> "BasesMask":
        """Parse a bases mask string into a BasesMask object."""
        reads = []
        str_parts = mask_str.split(",")
        regex = r"(?P<letter>[yniu])(?P<num>[0-9]*)"
        for part in str_parts:
            pieces = []
            for match in re.finditer(regex, part, flags=re.I):
                letter = match.group("letter").lower()
                num_str = match.group("num")
                num = 1 if len(num_str) == 0 else int(num_str)
                if pieces and pieces[-1][0] == letter:
                    # Collapse same-letter adjacent pieces
                    pieces[-1] = (pieces[-1][0], pieces[-1][1] + num)
                else:
                    pieces.append((letter, num))
            reads.append(pieces)
        return cls(reads=reads)

    def __str__(self) -> str:
        """Convert to bcl2fastq-style mask string (comma-separated, lowercase)."""

        def format_piece(letter: str, num: int) -> str:
            if num == 0:
                return ""
            elif num == 1:
                return letter
            return f"{letter}{num}"

        return ",".join(
            "".join(format_piece(*piece) for piece in read) for read in self.reads
        )

    def to_override_cycles(self) -> str:
        """
        Convert to bcl-convert OverrideCycles format.

        bcl-convert uses semicolons between reads and uppercase letters:
        - Y = use bases for read
        - I = use bases for index
        - U = use bases for UMI
        - N = skip bases

        Example: Y151;I8;I8;Y151
        """

        def format_piece(letter: str, num: int) -> str:
            if num == 0:
                return ""
            upper = letter.upper()
            if num == 1:
                return upper
            return f"{upper}{num}"

        return ";".join(
            "".join(format_piece(*piece) for piece in read) for read in self.reads
        )

    @property
    def num_index_reads(self) -> int:
        """Count the number of index reads in the mask."""
        return sum(1 for read in self.reads if any(piece[0] == "i" for piece in read))

    @property
    def index_lengths(self) -> tuple[int, int]:
        """
        Return the lengths of index1 and index2.
        Returns (0, 0) if no indexes, (len1, 0) if single index.
        """
        len1, len2 = 0, 0
        index_num = 0
        for read in self.reads:
            is_index = any(piece[0] == "i" for piece in read)
            if is_index:
                read_len = sum(piece[1] for piece in read)
                index_num += 1
                if index_num == 1:
                    len1 = read_len
                elif index_num == 2:
                    len2 = read_len
        return (len1, len2)

    def adjust_for_barcode_lengths(self, bc1_len: int, bc2_len: int) -> "BasesMask":
        """
        Create a new BasesMask adjusted for actual barcode lengths.

        If barcodes are shorter than the index read, pads with 'n'.
        If barcodes are longer than the index read, truncates.
        """
        new_reads = []
        index_num = 0

        for read in self.reads:
            read_len = sum(piece[1] for piece in read)
            is_index = any(piece[0] == "i" for piece in read)

            if is_index:
                if any(piece[0] == "y" for piece in read):
                    raise ValueError(
                        f"Mixed read/index in barcode mask '{self}', "
                        "don't know how to deal with this"
                    )

                index_num += 1
                bc_len = bc1_len if index_num == 1 else bc2_len

                if bc_len >= read_len:
                    # Barcode fills or exceeds the read
                    new_reads.append([("i", read_len)])
                else:
                    # Barcode is shorter, pad with 'n'
                    new_reads.append([("i", bc_len), ("n", read_len - bc_len)])
            else:
                new_reads.append(read)

        return BasesMask(reads=new_reads)


def parser_setup():
    parser = argparse.ArgumentParser(
        description="Generate bcl-convert compatible sample sheets from processing.json"
    )
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
    parser.add_argument(
        "--mismatches",
        help="Barcode mismatches allowed, as 'index1,index2' (default: 1,1)",
    )
    parser.set_defaults(**SCRIPT_OPTIONS)
    return parser


def get_barcode_assignments(
    data: dict, reverse_barcode1: bool, reverse_barcode2: bool
) -> "list[dict]":
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


def make_samplesheet_header(
    name: str, date: str, mask: BasesMask, mismatches: str
) -> str:
    """
    Generate bcl-convert v2 format sample sheet header.

    Args:
        name: Project/investigator name
        date: Date string
        mask: BasesMask object for the override cycles
        mismatches: Comma-separated mismatches for index1,index2 (e.g., "1,1")
    """
    mismatch_parts = mismatches.split(",")
    mismatch1 = mismatch_parts[0] if len(mismatch_parts) > 0 else "1"
    mismatch2 = mismatch_parts[1] if len(mismatch_parts) > 1 else mismatch1

    # Build mismatch settings - only include Index2 if there are 2 index reads
    num_indexes = mask.num_index_reads
    if num_indexes >= 2:
        mismatch_settings = (
            f"BarcodeMismatchesIndex1,{mismatch1}\nBarcodeMismatchesIndex2,{mismatch2}"
        )
    elif num_indexes == 1:
        mismatch_settings = f"BarcodeMismatchesIndex1,{mismatch1}"
    else:
        mismatch_settings = ""

    template = textwrap.dedent("""\
    [Header]
    Investigator Name,{name}
    Project Name,{name}
    Experiment Name,{name}
    Date,{date}
    Workflow,GenerateFASTQ

    [Settings]
    OverrideCycles,{override_cycles}
    {mismatch_settings}

    [Data]
    Lane,Sample_ID,index,index2
    """)
    return template.format(
        date=date,
        name=name,
        override_cycles=mask.to_override_cycles(),
        mismatch_settings=mismatch_settings,
    )


def group_assignments(assignments: "list[dict]") -> "dict[tuple[int, int], list[dict]]":
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


def write_samplesheets(
    name: str,
    filename_template: str,
    date: str,
    root_mask: str,
    assignments: list[dict],
    mismatches: str = "1,1",
) -> None:
    """Write out the sample sheets in bcl-convert v2 format."""
    mask = BasesMask.parse(root_mask)
    max_bclen1, max_bclen2 = mask.index_lengths

    # Trim barcodes to make sure they fit in the read
    for assign in assignments:
        assign["barcode1"] = assign["barcode1"][:max_bclen1]
        assign["barcode2"] = assign["barcode2"][:max_bclen2]

    groups = group_assignments(assignments)

    for barcode_lengths, assigns in groups.items():
        adjusted_mask = mask.adjust_for_barcode_lengths(*barcode_lengths)
        header = make_samplesheet_header(name, date, adjusted_mask, mismatches)
        body = make_samplesheet_body(assigns)
        samplesheet_contents = header + body
        filename = filename_template.format(mask=str(adjusted_mask))
        print(
            f"Writing {filename} with OverrideCycles={adjusted_mask.to_override_cycles()}"
        )
        with open(filename, "w") as f:
            f.write(samplesheet_contents)


def make_samplesheet_body(barcode_assignments: "list[dict]") -> str:
    """Create samplesheet text from assignments"""
    lines = []
    for ba in barcode_assignments:
        # bcl-convert format: Lane,Sample_ID,index,index2
        line = ",".join(
            [
                str(ba["lane"]),
                ba["sample"],
                str(ba["barcode1"]),
                str(ba["barcode2"]),
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
        mismatches=poptions.mismatches,
    )


if __name__ == "__main__":
    main()
