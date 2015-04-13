#!/bin/env python

import sys
import pysam
import argparse

parser = argparse.ArgumentParser(description="Set read pair status")
parser.add_argument('infile', metavar='IN_BAM', type=str,
                    help="The BAM file to read from. If '-', read from STDIN")
parser.add_argument('outfile', metavar='OUT_BAM', type=str,
                    help="The BAM file to write to. If '-', write to STDOUT")
parser.add_argument('--umi', action="store_true",
                    help="Extract and set UMI adapter information")

poptions = parser.parse_args()

unfiltered_reads = pysam.Samfile(poptions.infile, "rb")
filtered_reads = pysam.Samfile(poptions.outfile, "wbu", template = unfiltered_reads)

while(1):

  # pull the reads

  try:

    read1 = unfiltered_reads.next()
    read2 = unfiltered_reads.next()
    (read1, read2) = (read1, read2) if read1.is_read1 else (read2, read1)

  except:

    break

  # strip off the umi, and place it in a custom tag (if it exists)

  if poptions.umi:
      try:

        read1_umi_loc = read1.qname.index('#')
        read2_umi_loc = read2.qname.index('#')

      except:

        pass

      else:

        read1.setTag("XD", read1.qname[read1_umi_loc+1:])
        read1.qname =      read1.qname[:read1_umi_loc]

        read2.setTag("XD", read2.qname[read2_umi_loc+1:])
        read2.qname =      read2.qname[:read2_umi_loc]

  # filtering code

  try:

    # must be mapped to opposite strands (F-R conformation)
    if read1.is_reverse == read2.is_reverse: raise

    # both reads must have mapq greater than 0
    if read1.mapq == 0 or read2.mapq == 0: raise

    # each read must map to the same contig
    if read1.tid != read2.tid: raise

    # each read pair must be have an insert length of <=750 nt
    if read1.isize > 750 or read2.isize > 750: raise

  except:

    # failed a test above, not properly paired
    read1.flag &= ~(1<<1)
    read2.flag &= ~(1<<1)

  else:

    # pass all criteria, properly paired
    read1.flag |= (1<<1)
    read2.flag |= (1<<1)

  finally:

    # write to file
    filtered_reads.write(read1)
    filtered_reads.write(read2)

# clean-up and close files
unfiltered_reads.close()
filtered_reads.close()

