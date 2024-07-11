# This is a quick script to split up FASTQ files by barcode given
# Used to rescue tags from undeterminde state

import argparse
import errno
import itertools
import json
import logging
import operator
import os
import re
import subprocess
import sys

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

lengths = set([])


def parseArgs():
    parser = argparse.ArgumentParser(description="Split up fastq files by barcode")
    parser.add_argument(
        "--mismatches", type=int, default=0, help="number of mismatches"
    )
    parser.add_argument(
        "--processing",
        dest="processing_file",
        required=True,
        help="processing.json to use (mandatory)",
    )
    parser.add_argument(
        "--suffix", dest="suffix", default="", help="suffix to add to sample names"
    )
    parser.add_argument(
        "--lane", dest="lane", type=int, default=1, help="Lane to process (default 1)"
    )
    parser.add_argument(
        "--autosuffix",
        action="store_true",
        default=False,
        help="Automatically guess a suffix name",
    )
    parser.add_argument("--outdir", dest="outdir", default=".", help="Output directory")
    parser.add_argument(
        "--dry-run",
        dest="dry_run",
        default=False,
        action="store_true",
        help="Do not actually demultiplex",
    )
    parser.add_argument(
        "--ignore_failed_lanes",
        default=False,
        action="store_true",
        help="Ignore any lanes marked as failed in processing",
    )
    parser.add_argument(
        "--debug", dest="debug", action="store_true", help="Output debugging info"
    )
    parser.add_argument("infile", nargs="+")

    args = parser.parse_args()
    return args


def mismatch(word, mismatches):
    """Generator for mismatches
    returns original string + whatever variations have at most
    [mismatches] hamming distance"""
    for d in range(mismatches + 1):
        for locs in itertools.combinations(range(len(word)), d):
            this_word = [[char] for char in word]
            for loc in locs:
                orig_char = word[loc]
                this_word[loc] = [base for base in "ACGTN" if base != orig_char]
            for poss in itertools.product(*this_word):
                yield "".join(poss)


def guess_suffix(filename):
    file = os.path.basename(filename)

    nextseq_format = re.compile(
        r"""
    ^                       # Start
    Undetermined            #
    _ S0                    # index
    _ L (?P<lane> \d{3} )   # lane number
                            # Start suffix
    _ R (?P<read>   \d    ) # Read number
    _   (?P<count>  \d{3} ) # Count
                            # End suffix
    .fastq.gz
    $
    """,
        re.X,
    )

    match = nextseq_format.search(file)
    if match:
        suffix = "_R%s_%s" % (match.group("read"), match.group("lane"))
        return suffix

    hiseq_format = re.compile(
        r"""
    ^                  # Start
    lane \d+ _         # lane number
    Undetermined       #
    _ L (\d{3})        # lane number again
    (?P<suffix>        # Start suffix
      _ R \d           # Read number
      _   \d{3}        # Count
    )                  # End suffix
    .fastq.gz
    $
    """,
        re.X,
    )

    match = hiseq_format.search(file)
    if match:
        return match.group("suffix")

    return ""


def parse_processing_file(
    file, mismatches, suffix, lane, outdir, ignore_failed_lanes=False
):
    """Decode processing file into barcode->file assignments"""
    barcodes = {}
    labels = {}
    with open(file) as data_file:
        data = json.load(data_file)

    run_type = data["flowcell"]["run_type"]
    index_len = data["flowcell"]["index_length"]
    # Only some flowcell types need to treat different lanes differently
    if run_type == "NextSeq 500":
        lane_libraries = data["libraries"]
    elif run_type == "HISEQ V4":
        lane_libraries = [lib for lib in data["libraries"] if lib["lane"] == lane]
    elif run_type == "HiSeq 4000":
        lane_libraries = [lib for lib in data["libraries"] if lib["lane"] == lane]
    # TODO: Is this always correct?
    elif run_type.startswith("Novaseq 6000"):
        lane_libraries = [lib for lib in data["libraries"] if lib["lane"] == lane]
    else:
        logging.warn("Run type %s not supported; using all libraries" % run_type)
        lane_libraries = data["libraries"]

    for library in lane_libraries:
        if library.get("alignments", []):
            label = library["alignments"][0]["sample_name"]
        else:
            label = "%s_%s_L%03d" % (
                library["samplesheet_name"],
                library["barcode_index"],
                library["lane"],
            )

        if ignore_failed_lanes and library["failed"]:
            logging.info("Ignoring failed library %s" % label)
            continue

        project_dir = "Project_%s" % library["project"]
        sample_dir = "Sample_%s" % library["samplesheet_name"]
        library_dir = os.path.join(outdir, project_dir, sample_dir)
        outfile_name = os.path.join(library_dir, "%s%s.fastq.gz" % (label, suffix))

        try:
            os.makedirs(library_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        barcode_indices = library["barcode_index"].split("-")
        barcode1 = barcode_indices[0]
        barcode2 = barcode_indices[1] if len(barcode_indices) > 1 else ""

        # If any barcodes are longer than the index length
        # Trim those barcodes down to match.
        if len(barcode1) > index_len:
            barcode1 = barcode1[:index_len]
        if len(barcode2) > index_len:
            barcode2 = barcode2[:index_len]

        lengths.add((len(barcode1), len(barcode2)))

        for b1 in mismatch(barcode1, mismatches):
            for b2 in mismatch(barcode2, mismatches):
                barcode = (b1, b2)
                # TODO: This can be smarter
                if barcode in barcodes:
                    logging.error(
                        "Barcode %s already taken, lower --mismatches! (taken by %s+%s)"
                        % (barcode, barcode1, barcode2)
                    )
                    sys.exit(1)
                barcodes[barcode] = label

        labels[label] = {"filtered": 0, "unfiltered": 0, "total": 0}
        # TODO: Warning! this will overwrite files!
        outfile = open(outfile_name, "wb")
        labels[label]["fh"] = outfile
        labels[label]["out"] = subprocess.Popen(
            ["gzip", "-7"], stdout=outfile, stdin=subprocess.PIPE
        )

    logging.info(
        "Mapping %d barcodes to %s libraries" % (len(barcodes), len(lane_libraries))
    )
    logging.debug(barcodes)

    return barcodes, labels


def split_file(filename, barcodes, labels):
    """Demultiplex a file"""
    barcode_re = re.compile(
        r"""
            [012]:                  #
            ([YN]):                 # Fail/pass chastity filtering
            [01]:                   #
            ( [AGCTN] {6,20} )       # First barcode
            \+?                     # Optional separator (+)
            ( [AGCTN] {6,20} )?      # Optionally, second barcode
            $
            """,
        re.X,
    )

    tally = 0
    logging.info("Demultiplexing file: %s" % filename)

    if filename.endswith(".gz"):
        parsein = subprocess.Popen(
            ["zcat", filename], stdout=subprocess.PIPE, universal_newlines=True
        )
    else:
        parsein = subprocess.Popen(
            ["cat", filename], stdout=subprocess.PIPE, universal_newlines=True
        )

    for record, seq, qual in FastqGeneralIterator(parsein.stdout):
        tally += 1
        match = barcode_re.search(record)

        if not match:
            logging.error("Could not match %s" % record)
            logging.error(str(seq))
            logging.error("Record %d in %s" % (tally, filename))
            sys.exit(1)

        matches = match.groups()
        filter_str = matches[0]
        barcode1 = matches[1]
        barcode2 = ""
        label = None
        if len(matches) > 2:
            barcode2 = matches[2]
            if barcode2 is None:
                barcode2 = ""
        matched = False
        for fmt in lengths:
            barcode = (barcode1[: fmt[0]], barcode2[: fmt[1]])
            if barcode in barcodes:
                label = barcodes[barcode]
                matched = True
                break

        if matched:
            labels[label]["total"] += 1

            # Replace recorded barcode
            sep_index = record.rfind(":")
            record = record[: sep_index + 1] + barcode1 + "+" + barcode2
            # write to FASTQ
            text = bytes("@" + record + "\n" + seq + "\n+\n" + qual + "\n", "UTF-8")
            labels[label]["out"].stdin.write(text)

            if filter_str == "Y":
                labels[label]["filtered"] += 1
            else:
                labels[label]["unfiltered"] += 1

    # Wait for all subprocesses to finish
    parsein.communicate()


def main(argv):
    args = parseArgs()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format=log_format)
    else:
        logging.basicConfig(level=logging.INFO, format=log_format)

    global lengths
    lengths = set([])

    logging.info("File(s): %s" % args.infile)
    logging.info("OutDir:  %s" % args.outdir)
    logging.info("JSON:    %s" % args.processing_file)

    if args.autosuffix:
        args.suffix = guess_suffix(args.infile[0])
        logging.info("--autosuffix, guessing suffix as %s" % args.suffix)

    barcodes, labels = parse_processing_file(
        args.processing_file,
        args.mismatches,
        args.suffix,
        args.lane,
        args.outdir,
        ignore_failed_lanes=args.ignore_failed_lanes,
    )

    if args.dry_run:
        logging.info("Dry run, exiting")
        sys.exit(0)

    for filename in args.infile:
        split_file(filename, barcodes, labels)

    print("Barcode matching tallies:")

    for label, info in labels.items():
        print("%s\t%s" % (label, str(info)))

    for label, info in labels.items():
        info["out"].communicate()
        info["fh"].close()


if __name__ == "__main__":
    main(sys.argv)
