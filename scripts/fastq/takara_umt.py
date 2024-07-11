#!/usr/bin/env python3
import argparse
import logging

from Bio import SeqIO

LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

UMI_LEN = 8  # Fixed length for takara RNA UMTs
STEM_LEN = 6
TRIM_LEN = UMI_LEN + STEM_LEN


def parse_args():
    """Just parse the args"""
    parser = argparse.ArgumentParser(
        description="Annotate read names with Takara v3 UMT"
    )
    parser.add_argument(
        "--readlength",
        required=True,
        type=int,
        help="The length of each fastq file - used for trimming R1",
    )
    parser.add_argument("r1_fastq")
    parser.add_argument("r2_fastq")
    parser.add_argument("out_r1")
    parser.add_argument("out_r2")

    args = parser.parse_args()
    return args


def attach_umt(r1, r2, maxlen):
    """Attach UMT to r1 & r2 and remove from sequences"""
    # Put UMT in names
    umt = r2.seq[:UMI_LEN]
    umt_add = "#%s" % (umt)
    r1.id += umt_add
    r1.name = ""
    r1.description = " ".join(r1.description.split()[1:])
    r2.id += umt_add
    r2.description = " ".join(r2.description.split()[1:])
    r2.name = ""

    # Normal sequencing adapters have already been trimmed before this program
    # If adapter trimming has been performed, R1 ends in the UMT part
    # We can detect this by checking the length of the read against maxlen
    r1_trimmed = False
    if len(r1.seq) < maxlen:
        r1_trimmed = True
        r1 = r1[:-TRIM_LEN]

    # Remove from r2
    r2 = r2[TRIM_LEN:]
    return (r1, r2, r1_trimmed)


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)

    fragment_count = 0
    trim_count = 0
    removed_count = 0
    with open(args.r1_fastq) as r1_in, open(args.r2_fastq) as r2_in, open(
        args.out_r1, "wt"
    ) as r1_out, open(args.out_r2, "wt") as r2_out:
        r1_seq_io = SeqIO.parse(r1_in, "fastq")
        r2_seq_io = SeqIO.parse(r2_in, "fastq")

        for r1, r2 in zip(r1_seq_io, r2_seq_io):
            fragment_count += 1
            (r1, r2, trimmed) = attach_umt(r1, r2, args.readlength)
            if trimmed:
                trim_count += 1
            if len(r1.seq) == 0 or len(r2.seq) == 0:
                removed_count += 1
                continue
            r1_out.write(r1.format("fastq"))
            r2_out.write(r2.format("fastq"))

        logging.info("Run complete.")
        logging.info("Total Fragments: %d", fragment_count)
        logging.info("Trimming performed on: %d", trim_count)
        logging.info("Fragments omitted due to 0-length: %d", removed_count)


if __name__ == "__main__":
    main()
