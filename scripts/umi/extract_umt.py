#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import itertools
import logging
import string
import sys

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

UMI_LEN = 6

stems = ["TCAGTAGCTCA", "CAGTAGCTCA", "AGTAGCTCA", "GTAGCTCA"]
stem_lengths = [11, 10, 9, 8]

mismatched_stems = set()


def parseArgs():
    parser = argparse.ArgumentParser(description="Annotate read names with UMT")
    parser.add_argument(
        "--mismatches", type=int, default=1, help="number of mismatches"
    )
    parser.add_argument("r1_fastq")
    parser.add_argument("r2_fastq")
    parser.add_argument("out_r1")
    parser.add_argument("out_r2")

    args = parser.parse_args()
    return args


# Generator for mismatches
# returns original string + whatever variations have at most [mismatches]
# hamming distance
def mismatch(word, mismatches):
    for d in range(mismatches + 1):
        for locs in itertools.combinations(range(len(word)), d):
            thisWord = [[char] for char in word]
            for loc in locs:
                origChar = word[loc]
                thisWord[loc] = [l for l in "ACGTN" if l != origChar]
            for poss in itertools.product(*thisWord):
                yield "".join(poss)


# Finds the length of a stem in a sequence
# Relies on global variables:
# mismatched_stems, stem_lengths, and UMI_LEN
def find_stem_len(read):
    for len in stem_lengths:
        if str(read.seq[UMI_LEN : len + UMI_LEN]) in mismatched_stems:
            return len
    return 0


def attach_umt(r1, r2):
    r1_len = find_stem_len(r1)
    if not r1_len:
        return (None, None)

    r2_len = find_stem_len(r2)
    if not r2_len:
        return (None, None)

    # Put UMT in names
    umt_add = "#%s:%s" % (r1.seq[:UMI_LEN], r2.seq[:UMI_LEN])

    r1.id += umt_add
    r1.name = ""
    r1.description = " ".join(r1.description.split()[1:])

    r2.id += umt_add
    r2.description = " ".join(r2.description.split()[1:])
    r2.name = ""

    # Check for presence of UMT in mate - this indicates a short fragment that needs trimmed

    # Save stem & UMT for trimming use
    stem1 = r1[: UMI_LEN + r1_len]
    stem2 = r2[: UMI_LEN + r2_len]

    # Trim UMT & stem out of start of read
    r1 = r1[UMI_LEN + r1_len :]
    r2 = r2[UMI_LEN + r2_len :]

    # Trim ends, if necessary
    rev_stem1 = str(stem1.seq.reverse_complement())
    rev_stem2 = str(stem2.seq.reverse_complement())
    str_r1 = str(r1.seq)
    str_r2 = str(r2.seq)
    for x in range(UMI_LEN, -12, -1):
        x1 = x + r2_len
        x2 = x + r1_len
        if str_r1[-x1:] == rev_stem2[:x1] and str_r2[-x2:] == rev_stem1[:x2]:
            r1 = r1[:-x1]
            r2 = r2[:-x2]
            return (r1, r2)
    # Check specifically for "we didn't trim off one base of adapter sequence"
    x1 = UMI_LEN + r2_len
    x2 = UMI_LEN + r1_len
    if (
        str_r1[-x1 - 1 : -1] == rev_stem2[:x1]
        and str_r2[-x2 - 1 : -1] == rev_stem1[:x2]
    ):
        r1 = r1[: -x1 - 1]
        r2 = r2[: -x2 - 1]

    return (r1, r2)


# Modifies the global var mismatched_stems
def setup_mismatches(num_mismatches):
    for stem in stems:
        for mm in mismatch(stem, num_mismatches):
            mismatched_stems.add(mm)


def main(argv):
    args = parseArgs()
    logging.basicConfig(level=logging.WARN, format=log_format)
    setup_mismatches(args.mismatches)

    with open(args.r1_fastq) as r1_in, open(args.r2_fastq) as r2_in, open(
        args.out_r1, "wt"
    ) as r1_out, open(args.out_r2, "wt") as r2_out:
        r1_seqIO = SeqIO.parse(r1_in, "fastq")
        r2_seqIO = SeqIO.parse(r2_in, "fastq")
        try:
            while True:
                (r1, r2) = attach_umt(next(r1_seqIO), next(r2_seqIO))
                # Only write Fastq records for which we find stems
                if r1 is not None and r2 is not None:
                    r1_out.write(r1.format("fastq"))
                    r2_out.write(r2.format("fastq"))
        except StopIteration:
            logging.info("EOF reached")


if __name__ == "__main__":
    main(sys.argv)
