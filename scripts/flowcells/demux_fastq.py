# This is a quick script to split up FASTQ files by barcode given
# Used to rescue tags from undeterminde state

import sys, os, gzip, re, operator
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import logging

import argparse

import errno

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

def parseArgs():
    parser = argparse.ArgumentParser(description='Split up fastq files by barcode')
    parser.add_argument('--mismatches', type=int, default=0, help='number of mismatches')
    parser.add_argument('--processing', dest='processing_file', action='store', required=True,
            help='processing.json to use (mandatory)')
    parser.add_argument('--suffix', dest='suffix', default='', help='suffix to add to sample names')
    parser.add_argument('--lane', dest='lane', type=int, default=1, help='Lane to process (default 1)')
    parser.add_argument('--autosuffix', action="store_true", default=False, help='Automatically guess a suffix name')
    parser.add_argument('--outdir', dest='outdir', default='.', help='Output directory')
    parser.add_argument('--dry-run', dest='dry_run', default=False, action='store_true', help='Do not actually demultiplex')
    parser.add_argument('infile', nargs='+')

    args = parser.parse_args()
    return args



import itertools
# Generator for mismatches
# returns original string + whatever variations have at most [mismatches] hamming distance
def mismatch(word, mismatches):
    for d in range(mismatches+1):
        for locs in itertools.combinations(range(len(word)), d):
            thisWord = [[char] for char in word]
            for loc in locs:
                origChar = word[loc]
                thisWord[loc] = [l for l in "ACGTN" if l != origChar]
            for poss in itertools.product(*thisWord):
                yield "".join(poss)

def guess_suffix(filename):
    file = os.path.basename(filename)

    nextseq_format = re.compile(r"""
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
    """, re.X)

    match = nextseq_format.search(file)
    if match:
        suffix = "_R%s_%s" % ( match.group('read'), match.group('lane'))
        return suffix

    hiseq_format = re.compile(r"""
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
    """, re.X)

    match = hiseq_format.search(file)
    if match:
        return match.group('suffix')

    return ""

import json
def parse_processing_file(file, mismatches, suffix, lane, outdir):
    barcodes = dict()
    labels = dict()
    with open(file) as data_file:
        data = json.load(data_file)

    run_type = data['flowcell']['run_type']
    # Only some flowcell types need to treat different lanes differently
    if run_type == "NextSeq 500":
        lane_libraries = data['libraries']
    elif run_type == "HISEQ V4":
        lane_libraries = [ l for l in data['libraries'] if l['lane'] == lane ]
    else:
        logging.warn("Run type %s not supported; using all libraries" % run_type)
        lane_libraries = data['libraries']

    for library in lane_libraries:
        label = library['alignments'][0]['sample_name']
        project_dir = "Project_%s" % library['project']
        sample_dir = "Sample_%s" % library['samplesheet_name']
        library_dir = os.path.join(outdir, project_dir, sample_dir)
        outfile_name = os.path.join(library_dir, "%s%s.fastq.gz" % (label, suffix))

        try:
            os.makedirs(library_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


        barcode_indices = library['barcode_index'].split("-")
        barcode1 = barcode_indices[0]
        barcode2 = barcode_indices[1] if len(barcode_indices) > 1 else "" 

        lengths.add(( len(barcode1), len(barcode2) ))

        for b1 in mismatch(barcode1, mismatches):
            for b2 in mismatch(barcode2, mismatches):
                barcode = (b1, b2)
                # TODO: This can be smarter
                if barcode in barcodes:
                    logging.error("Barcode %s already taken, lower --mismatches! (taken by %s+%s)" % (barcode, barcode1, barcode2))
                    sys.exit(1)
                barcodes[barcode] = label

        labels[label] = { "filtered": 0, "unfiltered": 0, "total": 0 }
        # TODO: Warning! this will overwrite files!
        labels[label]["outfile"] = gzip.open(outfile_name, 'wb')

    logging.info("Mapping %d barcodes to %s libraries" % (len(barcodes), len(lane_libraries)))
    logging.debug(barcodes)
        
    return barcodes, labels

def split_file(filename, barcodes, labels):

    barcode_re = re.compile(r"""
            [012]:                  #
            ([YN]):                 # Fail/pass chastity filtering
            [01]:                   #
            ( [AGCTN] {6,8} )       # First barcode
            \+?                     # Optional separator (+)
            ( [AGCTN] {6,8} )?      # Optionally, second barcode
            $
            """, re.X)

    tally = 0
    logging.info("Demultiplexing file: %s" % filename)

    if filename.endswith('.gz'):
        parsein = gzip.open(filename)
    else:
        parsein = open(filename, 'rU')

    for record, seq, qual in FastqGeneralIterator(parsein):
        tally += 1
        match = barcode_re.search(record)

        if not match:
            debug.error("Could not match %s" % record)
            debug.error(str(seq))
            debug.error("Record %d in %s" % (tally, filename))
            sys.exit(1)

        matches = match.groups()
        filter = matches[0]
        barcode1 = matches[1]
        if (len(matches) > 2):
            barcode2 = matches[2]
        if barcode2 == None:
            barcode2 = ""
        matched=False
        for format in lengths:
            barcode = ( barcode1[:format[0]], barcode2[:format[1]] )
            if barcode in barcodes:
                label = barcodes[barcode]
                matched=True
                break

        if matched:
            labels[label]['total'] += 1

            # Replace recorded barcode
            sepIndex = record.rfind(':')
            record = record[:sepIndex + 1] + barcode1 + "+" + barcode2
            #write to FASTQ
            labels[label]['outfile'].write('@%s\n%s\n+\n%s\n' % (record, seq, qual))

            if filter == "Y":
                labels[label]["filtered"] += 1
            else:
                labels[label]['unfiltered'] += 1

    parsein.close()


def main(argv):
    args = parseArgs()

    logging.basicConfig(level=logging.INFO, format=log_format)

    global lengths
    lengths = set([])

    logging.info("File(s): %s" % args.infile)
    logging.info("OutDir:  %s" % args.outdir)
    logging.info("JSON:    %s" % args.processing_file)

    if args.autosuffix:
        args.suffix = guess_suffix(args.infile[0])
        logging.info("--autosuffix, guessing suffix as %s" % args.suffix)

    barcodes, labels = parse_processing_file(args.processing_file, args.mismatches, args.suffix, args.lane, args.outdir)

    if args.dry_run:
        logging.info("Dry run, exiting")
        sys.exit(0)


    for filename in args.infile:
        split_file(filename, barcodes, labels) 

    print("Barcode matching tallies:")

    for label, info in labels.iteritems():
        print "%s\t%s" % (label, str(info))

    for label, info in labels.iteritems():
        info['outfile'].close()

if __name__ == "__main__":
    main(sys.argv)
