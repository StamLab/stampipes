#!/usr/bin/env python3
"""
filter_reads.py - Set SAM flag 0x200 (QC-fail) for reads failing various criteria.

Criteria:
 * Input file must be non-marked.
 * Read must be end-to-end mapped without soft clipping or indels and meet cutoffs for MAPQ and NM
 * Separately flags reads that do not map to the nuclear genome as 0x1000
 * PE reads must have both reads present and meeting all above criteria. Additionally, they must
   face each other on the same reference and have an insert length wit
"""

import argparse
import logging
import sys

import pysam

"""
Exception when a bad read is found
"""


class ReadError(Exception):
    pass


"""
Looks for the UMI embeded in the read name, places it in a tag and trims the read name
"""


def parse_umi(read):
    try:
        umi_loc = read.query_name.index("#")
    except Exception:
        pass
    else:
        read.set_tag("XD", read.query_name[umi_loc + 1 :])
        read.query_name = read.query_name[:umi_loc]

    return read


"""
General function to set the flag field
"""


def set_read_flag(read, flag, mark):
    if mark:
        read.flag |= 1 << flag
    else:
        read.flag &= ~(1 << flag)

    return read


def set_proper_pair(read, mark=True):
    return set_read_flag(read, 1, mark)


def set_qc_fail(read, mark=True):
    return set_read_flag(read, 9, mark)


def set_nonnuclear(read, mark=True):
    return set_read_flag(read, 12, mark)


"""
General function to check whether a individual read passes QC
Checks mapq, mismatches
Throws a read exception if fails QC check
"""


def validate_read(read, min_mapq=1, max_mismatches=2):
    if read.mapping_quality < min_mapq:
        raise ReadError("Read MAPQ < %d" % min_mapq)
    if read.is_unmapped:
        raise ReadError("Read not mapped")
    if read.get_tag("NM") > max_mismatches:
        raise ReadError("Read mismatches > %d" % max_mismatches)

    return read


parser = argparse.ArgumentParser(
    prog="filter_reads",
    description="manual corrects the flags in a single- or pair-end BAM alignment file",
)
parser.add_argument(
    "raw_alignment", type=str, help="Inupt raw alignment file (must be sorted by name"
)
parser.add_argument(
    "filtered_alignment",
    type=str,
    help="Output filtered alignment file (sorted by name)",
)
parser.add_argument("nuclear_chr", type=str, help="List of nuclear chromosomes to use")
parser.add_argument(
    "--min_mapq",
    action="store",
    type=int,
    default=10,
    help="Reads must have at least this MAPQ to pass filter [%(default)s]",
)
parser.add_argument(
    "--max_mismatches",
    action="store",
    type=int,
    default=2,
    help="Maximum mismatches to pass filter [%(default)s]",
)
parser.add_argument(
    "--max_insert_size",
    action="store",
    type=int,
    default=750,
    help="Maximum insert size to pass filter [%(default)s]",
)
parser.add_argument(
    "--verbosity",
    action="store",
    type=int,
    default=50,
    help="Verbosity (50 = quiet, 0 = loud) [%(default)s]",
)
args = parser.parse_args()

logging.basicConfig(stream=sys.stdout, level=args.verbosity)

raw_alignment = pysam.AlignmentFile(args.raw_alignment, "rb")
filtered_alignment = pysam.AlignmentFile(
    args.filtered_alignment, "wbu", template=raw_alignment
)
nuclear_chrs = [line.rstrip("\n") for line in open(args.nuclear_chr)]

raw_reads = raw_alignment.fetch(until_eof=True)

read1 = None
read2 = None

qc_fail = False
proper_pair = False

while 1:
    try:
        if not read1:
            read1 = parse_umi(next(raw_reads))
    except Exception:
        break

    try:
        read2 = parse_umi(next(raw_reads))
    except Exception:
        read2 = None

    # Continue in pair-end mode if their is two reads that are paired and that they have the same name

    if (
        (read1 and read2)
        and (read1.is_paired and read2.is_paired)
        and (read1.query_name == read2.query_name)
    ):
        (read1, read2) = (read1, read2) if read1.is_read1 else (read2, read1)

        try:
            # Check if the individual reads pass QC; throws exception if not

            validate_read(read1, args.min_mapq, args.max_mismatches)
            validate_read(read2, args.min_mapq, args.max_mismatches)

            # Read pair must be in F-R configuration

            if read1.is_reverse == read2.is_reverse:
                raise ReadError("Mates cannot align to same strand!")

            # Both reads must map to same contig

            if read1.reference_id != read2.reference_id:
                raise ReadError("Mates must align to the same reference contig!")

            # Insert size must be greater than 0

            if read1.template_length == 0 or read2.template_length == 0:
                raise ReadError("Insert size cannot be 0!")

            # Insert sizes must add up to zero (one must be positive and the other negative)

            if read1.template_length + read2.template_length != 0:
                raise ReadError("Insert sizes must be equal!")

            # Insert sizes must less than the maximum

            if (
                abs(read1.template_length) > args.max_insert_size
                or read2.template_length > args.max_insert_size
            ):
                raise ReadError("Insert size > %d!" % args.max_insert_size)

        except ReadError as e:
            # If we get a read exception, then set
            # QC fail flag, and unset proper pair flag

            qc_fail = True
            proper_pair = False

            logging.debug(e)
            logging.debug(read1)
            logging.debug(read2)

        else:
            # No exception, then unset
            # QC fail flag, and unset proper pair flag

            qc_fail = False
            proper_pair = True

        finally:
            # Set the flags

            set_qc_fail(read1, qc_fail)
            set_proper_pair(read1, proper_pair)

            if read1.reference_id != -1 and read1.reference_name not in nuclear_chrs:
                set_nonnuclear(read1, True)

            set_qc_fail(read2, qc_fail)
            set_proper_pair(read2, proper_pair)
            if read2.reference_id != -1 and read2.reference_name not in nuclear_chrs:
                set_nonnuclear(read2, True)

            # Write file

            filtered_alignment.write(read1)
            filtered_alignment.write(read2)

            (read1, read2) = (None, None)

    # Failed pair-end test -- could be single-end

    else:
        try:
            validate_read(read1, args.min_mapq, args.max_mismatches)

            if read1.is_paired:
                if not read1.mate_is_unmapped:
                    raise ReadError("No mate found (incongruent flag)!")

                else:
                    raise ReadError("No mate found!")

        except ReadError as e:
            qc_fail = True

            logging.debug(e)
            logging.debug(read1)

        else:
            qc_fail = False

        finally:
            set_qc_fail(read1, qc_fail)
            set_proper_pair(read1, False)

            if read1.reference_id != -1 and read1.reference_name not in nuclear_chrs:
                set_nonnuclear(read1, True)

            filtered_alignment.write(read1)

            (read1, read2) = (read2, None)

# clean-up and close files
raw_alignment.close()
filtered_alignment.close()
