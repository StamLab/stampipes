#!/usr/bin/env python

import argparse
import sys
import pysam
from collections import defaultdict

UMI_TAG = "XD:Z"


def parser_setup():

    script_options = {
        "infile": "/dev/stdin",
        "outfile": "/dev/stdout",
        "histfile": None,
    }

    parser = argparse.ArgumentParser()

    parser.add_argument("--hist", dest="histfile",
                        help="Write histogram of duplicates to this file")

    parser.add_argument("-i", "--infile", dest="infile",
                        help="Read from this file")
    parser.add_argument("-o", "--outfile", dest="outfile",
                        help="Write to this file")

    parser.set_defaults(**script_options)
    parser.set_defaults(quiet=False, debug=False)

    return parser


# Takes a list of reads, returns one sorted such that the "best" read is at the
# top.  That means highest mapping quality, with ties broken by query_name (in
# lexicographic order)
def sortQuality(x):
    return (-1 * x.mapping_quality, x.query_name)


# Takes in a list of reads, returns a list of lists - one list for each unique
# read. These reads should have identical chromosome & reference-start fields.
# They are sorted so that the highest-mapping quality read is first in each
# sub-list.
def find_dups(reads):
    # Don't waste time if there's only one read at this position
    if len(reads) == 1:
        return [reads]

    lists = defaultdict(list)
    for r in reads:
        if r.has_tag(UMI_TAG):
            key = (r.template_length, r.get_tag(UMI_TAG))
        else:
            key = r.template_length
        lists[key].append(r)

    return [sorted(l, key=sortQuality) for l in lists.values()]


# Sets a read's duplicate flag, returns it
def set_dup(read, is_dup):
    if is_dup:
        read.flag |= 1024
    else:
        read.flag &= ~1024
    return read


class DupMarker():

    read_histo = defaultdict(int)
    pair_map = {}
    input = None
    output = None
    histo = None

    def __init__(self, input, output=None, histo=None):
        self.input = input
        self.output = output
        self.histo = histo

    def process_chunk(self, chunk):
        new_reads = []
        already_seen_reads = []
        for r in chunk:
            if r.query_name in self.pair_map:
                already_seen_reads.append(r)
            else:
                self.pair_map[r.query_name] = None
                new_reads.append(r)

        read_sets = find_dups(new_reads)

        for l in read_sets:
            if self.histo is not None:
                self.read_histo[len(l)] += 1
            # Mark all but the highest-quality read as duplicates
            l[0] = set_dup(l[0], False)
            self.pair_map[l[0].query_name] = False
            for r in l[1:]:
                r = set_dup(r, True)
                self.pair_map[r.query_name] = True

            # Print the read set
            if self.output is not None:
                for r in l:
                    self.output.write(r)

        # Print out already seen reads
        for r in already_seen_reads:
            r = set_dup(r, self.pair_map[r.query_name])
            del self.pair_map[r.query_name]
            if self.output is not None:
                self.output.write(r)

    def run(self):
        chr = None
        pos = None
        current_chunk = []
        for read in self.input:
            if read.reference_id != chr or read.reference_start != pos:
                self.process_chunk(current_chunk)
                current_chunk = []
                chr = read.reference_id
                pos = read.reference_start
            current_chunk.append(read)

        self.process_chunk(current_chunk)

        if self.histo is not None:
            for k in sorted(self.read_histo.keys()):
                self.histo.write("%s\t%s\n" % (k, self.read_histo[k]))


def main(args=sys.argv):
    parser = parser_setup()
    poptions = parser.parse_args()

    input = pysam.AlignmentFile(poptions.infile, 'r')
    output = pysam.AlignmentFile(poptions.outfile, 'wb0', template=input)
    histo = open(poptions.histfile, 'w') if poptions.histfile else None

    try:
        dupmarker = DupMarker(input=input, output=output, histo=histo)
        dupmarker.run()
    finally:
        input.close()
        output.close()
        histo.close()


# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell
# after importing without automatically running it
if __name__ == "__main__":
    main()
