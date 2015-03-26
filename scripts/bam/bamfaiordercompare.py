"""
This quick script uses pysam's capabilities to compare the ordering
of references in the header of a BAM file with the order in a FAI
file.  Prints "ORDERED" if they are equal and "UNORDERED" otherwise.
"""

import sys
import os

import pysam

def compare_bam_order(faifile, bamfile):
    fai = open(faifile, 'r').read()

    bam = pysam.Samfile( bamfile, "rb")
    bamorder = list(bam.references)
    faiorder = [line.split('\t')[0] for line in fai.split('\n') if line]

    if bamorder == faiorder:
        return True
    return False

def main(args = sys.argv):

    if len(sys.argv) < 3:
        print "USAGE: %s FAIFILE BAMFILE" % sys.argv[0]
        sys.exit(0)

    faifile = sys.argv[1]
    bamfile = sys.argv[2]

    if not os.path.exists(faifile):
        sys.stderr.write("FAIFILE %s does not exist\n" % faifile)
        sys.exit(1)

    if not os.path.exists(bamfile):
        sys.stderr.write("BAMFILE %s does not exist\n" % bamfile)
        sys.exit(1)

    ordered = compare_bam_order(faifile, bamfile)

    if ordered:
        sys.stdout.write("ORDERED\n")
    else:
        sys.stdout.write("UNORDERED\n")

if __name__ == "__main__":
    main()
