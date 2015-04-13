#!/usr/bin/env python

"""
Filters a BWA aligned BAM file to only include uniquely mapping, properly
paired reads.

Useful SAM flag reference: http://broadinstitute.github.io/picard/explain-flags.html
"""

import os, sys, logging, re
from pysam import Samfile
from optparse import OptionParser

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

usage = """usage:
%prog BAMFILE COUNTFILE [OPTIONS]"""

script_options = {
  "debug": False,
  "quiet": True,
}

def parser_setup():

    parser = OptionParser(usage=usage)

    parser.add_option("-q", "--quiet", dest="quiet", action="store_true", 
        help="Don't print info messages to standard out.")
    parser.add_option("-d", "--debug", dest="debug", action="store_true", 
        help="Print all debug messages to standard out.")

    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser

class BAMFilter(object):

    def __init__(self):

        self.max_mismatches = 2
        self.previous_read = None
        self.min_map_quality = 30

    def process_read(self, read, inbam): 

        tags = dict(read.tags)

        mapq = "mapq-%d" % read.mapq
        if not mapq in self.mapqcounts:
            self.mapqcounts[mapq] = 1
        else:
            self.mapqcounts[mapq] += 1

        samflag = "samflag-%d" % read.flag
        if not samflag in self.samflagcounts:
            self.samflagcounts[samflag] = 1
        else:
            self.samflagcounts[samflag] += 1

        readlength = "readlength-%d" % read.rlen
        if not readlength in self.readlengthcounts:
            self.readlengthcounts[readlength] = 1
        else:
            self.readlengthcounts[readlength] += 1

        # This might not be the most perfect indicator, but it will do for now
        # Must take place before minimum quality filter--most multiple matching
        # reads have mapq set to 0
        if "XT" in tags and tags["XT"] == "R":
            self.counts['mm'] += 1

        if read.is_qcfail:
            self.counts['qc-flagged'] += 1

        if read.is_unmapped:
            self.counts['nm'] += 1
            return False
        else:
            self.counts['all-aligned'] += 1

        # Figure out how many alignments aren't included because of mapq
        if self.min_map_quality > read.mapq:
            self.counts['all-mapq-filter'] += 1
            return

        # do not use reads with
        # 0x4 = read is unmapped
        # 0x8 = pair is unmapped
        if read.flag & 12:
            return False

        # only use reads with
        # 0x1 read paired
        # 0x2 read mapped in proper pair
        if not (read.flag & 1 and read.flag & 2):
            return False

        # Figure out how many alignments aren't included because of mapq
        if self.min_map_quality > read.mapq:
            self.counts['paired-mapq-filter'] += 1
            return

        # do not use reads with QC fail even if they pass all other checks
        # 0x512 QC Fail
        if read.flag & 512:
            self.counts['paired-aligned-qcfail'] += 1
            return

        self.counts['paired-aligned'] += 1
        self.counts['u'] += 1

        passreadlength = "aligned-readlength-%d" % read.rlen
        if not passreadlength in self.readlengthcounts:
            self.readlengthcounts[passreadlength] = 1
        else:
            self.readlengthcounts[passreadlength] += 1

        if read.is_qcfail:
            return False

        self.counts['u-pf'] += 1

        if "N" in read.seq:
            return False
        else:
            self.counts['u-pf-n'] += 1

        if tags['NM'] > self.max_mismatches:
            return False
        else:
            self.counts['u-pf-n-mm%d' % self.max_mismatches] += 1

        chr = inbam.getrname(read.rname)

        if not chr in self.chrcounts:
            self.chrcounts[chr] = 1
        else:
            self.chrcounts[chr] += 1

        if not "chrM" == chr:
            self.counts['u-pf-n-mm%d-mito' % self.max_mismatches] += 1

        self.previous_read = read

        return True

    def write_dict(self, countout, counts):

        for count in sorted(counts.keys()):
            countout.write("%s\t%d\n" % (count, counts[count]))

    def filter(self, infile, countfile):

        inbam = Samfile(infile, 'rb')

        count_labels = ['u', 'u-pf', 'u-pf-n', 'u-pf-n-mm%d' % self.max_mismatches,
          'u-pf-n-mm%d-mito' % self.max_mismatches, 'mm', 'nm', 'qc-flagged',
          'all-aligned', 'paired-aligned', 'all-mapq-filter', 'paired-mapq-filter', 'paired-aligned-qcfail']

        self.counts = dict([(label, 0) for label in count_labels])

        self.chrcounts = {}
        self.mapqcounts = {}
        self.samflagcounts = {}
        self.readlengthcounts = {}

        for read in inbam:
            self.process_read(read, inbam)

        countout = open(countfile, 'a')

        self.write_dict(countout, self.counts)
        self.write_dict(countout, self.chrcounts)
        self.write_dict(countout, self.mapqcounts)
        self.write_dict(countout, self.samflagcounts)
        self.write_dict(countout, self.readlengthcounts)

        countout.close()

def main(args = sys.argv):
    """This is the main body of the program that by default uses the arguments
from the command line."""

    parser = parser_setup()
    (poptions, pargs) = parser.parse_args()

    if poptions.quiet:
        logging.basicConfig(level=logging.WARNING, format=log_format)
    elif poptions.debug:
        logging.basicConfig(level=logging.DEBUG, format=log_format)
    else:
        # Set up the logging levels
        logging.basicConfig(level=logging.INFO, format=log_format)

    bamfile = pargs[0]
    # outfile = pargs[1]
    countfile = pargs[1]

    filter = BAMFilter()

    filter.filter(bamfile, countfile)

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()
