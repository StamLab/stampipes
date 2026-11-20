#!/usr/bin/env python3
"""
Generate bcl-convert compatible sample sheets from processing.json

This script creates sample sheets with the bcl-convert v2 format, including:
- OverrideCycles in [Settings] section (replaces bcl2fastq's --use-bases-mask)
- BarcodeMismatchesIndex1/2 in [Settings] section
- Proper [Data] section with Lane, Sample_ID, index, index2 columns

Barcode values come from a LIMS field policy (--barcode-source). That is
separate from RunInfo.xml IsReverseComplement, which controls how bcl-convert
interprets OverrideCycles for each physical index read (including padding
order). Do not conflate the two.

Usage:
  $0 -p processing.json --runinfo RunInfo.xml [--barcode-source reverse_sequence]
     [--mismatches 1,1] [--filename SampleSheet.withmask.{mask}.csv]
"""

import argparse
import datetime
import json
import re
import textwrap
import xml.etree.ElementTree as ElementTree
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

SCRIPT_OPTIONS = {
    "processing": "processing.json",
    "runinfo": "RunInfo.xml",
    "barcode_source": "sequence",
    "filename": "SampleSheet.withmask.{mask}.csv",
    "mismatches": "1,1",
}

BARCODE_SOURCES = ("sequence", "reverse_sequence")

# LIMS barcode-category umi_method values whose UMI is read immediately after
# the barcode in Index Read 1 (e.g. IDT xGen UDI-UMI: 8bp i7 + 9bp UMI).
# Value is the UMI length in cycles.
INDEX1_UMI_LENGTHS = {
    "xgen-udi-umi": 9,
}
# Fallback keyed on LIMS barcode-category name, for categories that don't have
# a umi_method set yet. Fragile (human-editable name); prefer umi_method.
INDEX1_UMI_CATEGORY_LENGTHS = {
    "IDT xGen UDI-UMI Adapters": 9,
}


@dataclass(frozen=True)
class IndexReadInfo:
    """Physical index-read metadata from RunInfo.xml, in read order."""

    num_cycles: int
    is_reverse_complement: bool


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
    def index_active_states(self) -> list[bool]:
        """
        Return active state for each physical index position.

        A physical index slot is identified by containing any 'i' piece (including
        zero-length ones left by adjust_for_barcode_lengths).  A slot is *active*
        when it has at least one 'i' piece with length > 0, meaning bcl-convert
        will actually demultiplex on it.

        Returns a list ordered by physical read position, e.g. [True, False] means
        Index1 is active and Index2 is masked.
        """
        states = []
        for read in self.reads:
            if any(piece[0] == "i" for piece in read):
                active = any(piece[0] == "i" and piece[1] > 0 for piece in read)
                states.append(active)
        return states

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

    def adjust_for_barcode_lengths(
        self,
        bc1_len: int,
        bc2_len: int,
        index_reads: "list[IndexReadInfo] | None" = None,
        umi1_len: int = 0,
    ) -> "BasesMask":
        """
        Create a new BasesMask adjusted for actual barcode lengths.

        If barcodes are shorter than the index read, pads with 'n'.
        If barcodes are longer than the index read, truncates.

        For index reads marked IsReverseComplement in RunInfo.xml, bcl-convert
        reverses the OverrideCycles token before interpreting cycles. Padding
        must therefore lead on those reads (N2I8) rather than trail (I8N2).

        If umi1_len is set, the umi1_len cycles after barcode1 in Index Read 1
        are marked as UMI (I8U9, or I8U9N2 if cycles remain).
        """
        if index_reads is None:
            raise ValueError(
                "index_reads from RunInfo.xml is required to place OverrideCycles "
                "padding correctly for reverse-complement index reads"
            )

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

                if index_num >= len(index_reads):
                    raise ValueError(
                        f"Bases mask '{self}' has more index reads than "
                        f"RunInfo.xml ({len(index_reads)})"
                    )

                index_info = index_reads[index_num]
                index_num += 1
                bc_len = bc1_len if index_num == 1 else bc2_len

                if index_num == 1 and umi1_len:
                    if index_info.is_reverse_complement:
                        raise ValueError(
                            "UMI on a reverse-complement index read is not supported"
                        )
                    pad = read_len - bc_len - umi1_len
                    if pad < 0:
                        raise ValueError(
                            f"Index read 1 has {read_len} cycles, too short for a "
                            f"{bc_len}bp barcode plus {umi1_len}bp UMI"
                        )
                    pieces = [("i", bc_len), ("u", umi1_len)]
                    if pad:
                        pieces.append(("n", pad))
                    new_reads.append(pieces)
                elif bc_len >= read_len:
                    # Barcode fills or exceeds the read
                    new_reads.append([("i", read_len)])
                else:
                    # Barcode is shorter (or absent), pad with 'n'.
                    # Keep ('i', 0) as a zero-length marker so callers can
                    # distinguish a masked physical index slot from a non-index read.
                    pad = read_len - bc_len
                    if index_info.is_reverse_complement:
                        new_reads.append([("n", pad), ("i", bc_len)])
                    else:
                        new_reads.append([("i", bc_len), ("n", pad)])
            else:
                new_reads.append(read)

        if index_num != len(index_reads):
            raise ValueError(
                f"Bases mask '{self}' has {index_num} index reads but "
                f"RunInfo.xml has {len(index_reads)}"
            )

        return BasesMask(reads=new_reads)


def parse_runinfo(path: "str | Path") -> list[IndexReadInfo]:
    """
    Parse physical index-read metadata from RunInfo.xml.

    Returns IndexReadInfo entries in RunInfo read order. Missing
    IsReverseComplement attributes are treated as forward (False).
    """
    runinfo_path = Path(path)
    if not runinfo_path.is_file():
        raise FileNotFoundError(f"RunInfo.xml not found: {runinfo_path}")

    try:
        root = ElementTree.parse(runinfo_path).getroot()
    except ElementTree.ParseError as exc:
        raise ValueError(f"Malformed RunInfo.xml ({runinfo_path}): {exc}") from exc

    reads = root.find("Run")
    if reads is not None:
        reads = reads.find("Reads")
    if reads is None:
        reads = root.find("Reads")
    if reads is None:
        raise ValueError(f"No <Reads> section in RunInfo.xml: {runinfo_path}")

    index_reads: list[IndexReadInfo] = []
    for read in reads.findall("Read"):
        is_indexed = read.attrib.get("IsIndexedRead", "N").upper() == "Y"
        if not is_indexed:
            continue
        try:
            num_cycles = int(read.attrib["NumCycles"])
        except (KeyError, ValueError) as exc:
            raise ValueError(
                f"Index read missing valid NumCycles in {runinfo_path}"
            ) from exc
        is_rc = read.attrib.get("IsReverseComplement", "N").upper() == "Y"
        index_reads.append(
            IndexReadInfo(num_cycles=num_cycles, is_reverse_complement=is_rc)
        )

    if not index_reads:
        raise ValueError(f"No indexed reads found in RunInfo.xml: {runinfo_path}")

    return index_reads


def validate_runinfo_against_mask(
    mask: BasesMask, index_reads: list[IndexReadInfo]
) -> None:
    """Ensure RunInfo index cycle counts match the bases mask."""
    physical_lengths: list[int] = []
    for read in mask.reads:
        if any(piece[0] == "i" for piece in read):
            physical_lengths.append(sum(piece[1] for piece in read))

    if len(physical_lengths) != len(index_reads):
        raise ValueError(
            f"Bases mask has {len(physical_lengths)} index reads but "
            f"RunInfo.xml has {len(index_reads)}"
        )

    for idx, (mask_len, info) in enumerate(zip(physical_lengths, index_reads), start=1):
        if mask_len != info.num_cycles:
            raise ValueError(
                f"Index read {idx} length mismatch: bases mask has {mask_len} "
                f"cycles but RunInfo.xml has {info.num_cycles}"
            )


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
        "--runinfo",
        dest="runinfo",
        help=(
            "Path to RunInfo.xml (default: RunInfo.xml). Used for physical "
            "index-read cycle counts and IsReverseComplement padding order."
        ),
    )
    parser.add_argument(
        "--barcode-source",
        dest="barcode_source",
        choices=BARCODE_SOURCES,
        help=(
            "Which LIMS barcode field to write into the sample sheet "
            f"({', '.join(BARCODE_SOURCES)}; default: sequence). This is a "
            "LIMS representation policy, not bcl-convert read orientation."
        ),
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


def get_barcode_assignments(data: dict, barcode_source: str) -> "list[dict]":
    """
    Select sample-sheet barcode strings from LIMS fields.

    barcode_source chooses sequence vs reverse_sequence. That choice is
    independent of RunInfo.xml IsReverseComplement / OverrideCycles padding.
    barcode_check.py compares against raw BCL barcodes and intentionally uses
    a different representation; do not unify the two.
    """
    if barcode_source not in BARCODE_SOURCES:
        raise ValueError(
            f"barcode_source must be one of {BARCODE_SOURCES}, got {barcode_source!r}"
        )

    assignments = []

    for libdata in data["libraries"]:
        assignment = {
            "lane": libdata.get("lane"),
            "sample": libdata.get("samplesheet_name"),
            "barcode1": "",
            "barcode2": "",
            "umi1": 0,
        }
        if assignment["sample"] == "None":
            assignment["sample"] = "LANE%d" % libdata["id"]
        if libdata.get("barcode1") is not None:
            assignment["barcode1"] = libdata["barcode1"][barcode_source]
            umi_method = libdata["barcode1"].get("umi_method")
            category = libdata["barcode1"].get("category")
            assignment["umi1"] = INDEX1_UMI_LENGTHS.get(
                umi_method, INDEX1_UMI_CATEGORY_LENGTHS.get(category, 0)
            )
        if libdata.get("barcode2") is not None:
            assignment["barcode2"] = libdata["barcode2"][barcode_source]

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

    # Emit BarcodeMismatchesIndexN only for physically active index reads.
    # bcl-convert uses physical read numbering: Index1 = first index slot in
    # RunInfo.xml, Index2 = second.  Emitting a mismatch setting for a slot
    # whose OverrideCycles are all-N causes bcl-convert to error.
    active = mask.index_active_states
    mismatch_lines = []
    if len(active) >= 1 and active[0]:
        mismatch_lines.append(f"BarcodeMismatchesIndex1,{mismatch1}")
    if len(active) >= 2 and active[1]:
        mismatch_lines.append(f"BarcodeMismatchesIndex2,{mismatch2}")
    mismatch_settings = "\n".join(mismatch_lines)

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


def group_assignments(
    assignments: "list[dict]",
) -> "dict[tuple[int, int, int], list[dict]]":
    """Groups the barcode assignments by (barcode1, barcode2, index1 UMI) length"""
    barcode_length_combinations = defaultdict(list)

    def get_len(d):
        return 0 if d is None else len(d)

    for assignment in assignments:
        key = (
            get_len(assignment["barcode1"]),
            get_len(assignment["barcode2"]),
            assignment.get("umi1", 0),
        )
        barcode_length_combinations[key].append(assignment)
    return barcode_length_combinations


def write_samplesheets(
    name: str,
    filename_template: str,
    date: str,
    root_mask: str,
    assignments: list[dict],
    index_reads: list[IndexReadInfo],
    mismatches: str = "1,1",
) -> None:
    """Write out the sample sheets in bcl-convert v2 format."""
    mask = BasesMask.parse(root_mask)
    validate_runinfo_against_mask(mask, index_reads)
    max_bclen1, max_bclen2 = mask.index_lengths

    # Trim barcodes to make sure they fit in the read
    for assign in assignments:
        assign["barcode1"] = assign["barcode1"][:max_bclen1]
        assign["barcode2"] = assign["barcode2"][:max_bclen2]

    groups = group_assignments(assignments)

    for (bc1_len, bc2_len, umi1_len), assigns in groups.items():
        adjusted_mask = mask.adjust_for_barcode_lengths(
            bc1_len, bc2_len, index_reads=index_reads, umi1_len=umi1_len
        )
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


def main(args=None):
    parser = parser_setup()
    poptions = parser.parse_args(args)

    with open(poptions.processing) as process_json:
        data = json.load(process_json)

    index_reads = parse_runinfo(poptions.runinfo)
    assignments = get_barcode_assignments(data, poptions.barcode_source)
    mask = data["alignment_group"]["bases_mask"]
    write_samplesheets(
        name="Altius",
        filename_template=poptions.filename,
        date=str(datetime.date.today()),
        root_mask=mask,
        assignments=assignments,
        index_reads=index_reads,
        mismatches=poptions.mismatches,
    )


if __name__ == "__main__":
    main()
