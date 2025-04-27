"""
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3959825/
For peak calling and footprinting, we adjusted the read start sites to represent the center of the transposon binding event.
Previous descriptions of the Tn5 transposase show that the transposon binds as a dimer and inserts two adapters separated by
9 bps (main text ref. 11). Therefore, all reads aligning to the + strand were offset by +4 bps, and all reads
aligning to the – strand were offset −5 bps.

python shiftbam.py < marked.bam > shifted.bam 2> shifted.log
"""

import pysam
import csv
import sys

def shift_bam():
    """
    Because of vagaries with pysam, this uses stdin and stdout:
    https://pysam.readthedocs.io/en/latest/usage.html#using-streams
    """
    infile = pysam.AlignmentFile("-", "rb")
    outfile = pysam.AlignmentFile("-", "wb", template=infile)

    for n, read in enumerate(infile):
        original_start = read.reference_start
        # make the changes in the reference start
        if read.is_reverse:
            read.reference_start -= 5
        else:
            read.reference_start += 4

        # we can't go below 0
        if read.reference_start < 0:
            print("Adjusted {} from {} to {} in {}, but setting to 0".format(
                read.query_name, original_start, read.reference_start, read.reference_name), file=sys.stderr)
            read.reference_start = 0
        # ended up not using inferred length because sometimes it is None?
        if read.reference_start + read.query_length > infile.lengths[read.reference_id]:
            new_location = infile.lengths[read.reference_id] - read.query_length
            print("Adjusted {} from {} to {} in {}, but setting to {}".format(
                read.query_name, original_start, read.reference_start, read.reference_name, new_location), file=sys.stderr)
            read.reference_start = new_location

        outfile.write(read)

if __name__ == "__main__":

    shift_bam()
