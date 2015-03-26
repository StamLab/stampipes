#!/bin/bash

if [ "$#" -lt "4" ] ; then
    echo "Usage:  $0  input_F.gz input_R.gz output_F.gz output_R.gz" 1>&2
    exit 1
fi
readsF=$1
readsR=$2
outputF=$3
outputR=$4

outputdir=`dirname $outputF`
mkdir -p $outputdir

adapterFA=$STAMPIPES/data/adapters/IlluminaAdapter_min40bp_revcomp.fa

trim-adapters-illumina -f $adapterFA --length=$READLENGTH $readsF $readsR $outputF $outputR


#  usage: fastq-mcf [options] <adapters.fa> <reads.fq> [mates1.fq ...] 
#  
#  Detects levels of adapter presence, computes likelihoods and
#  locations (start, end) of the adapters.   Removes the adapter
#  sequences from the fastq file, and sends it to stdout.
#  
#  Stats go to stderr, unless -o is specified.
#  
#  If you specify multiple 'paired-end' inputs, then a -o option is
#  required for each.  IE: -o read1.clip.q -o read2.clip.fq
#  
#  Options:
#  	-h      This help
#  	-o FIL  Output file (stats to stdout)
#  	-s N.N  Log scale for clip pct to threshold (2.5)
#  	-t N    % occurance threshold before clipping (0.25)
#  	-m N    Minimum clip length, overrides scaled auto (1)
#  	-p N    Maximum adapter difference percentage (20)
#  	-l N    Minimum remaining sequence length (15)
#  	-L N    Maximum sequence length (none)
#  	-k N    sKew percentage causing trimming (2)
#  	-q N    quality threshold causing trimming (10)
#  	-w N    window-size for quality trimming (1)
#  	-f      force output, even if not much will be done
#  	-0      Set all trimming parameters to zero
#  	-U|u    Force disable/enable illumina PF filtering
#  	-P N    phred-scale (64)
#  	-x N    'N' (Bad read) percentage causing trimming (10)
#  	-R      Don't remove N's from the fronts/ends of reads
#  	-n      Don't clip, just output what would be done
#  	-C N    Number of reads to use for subsampling (100k)
#  	-S FIL  Save clipped reads to file
#  	-d      Output lots of random debugging stuff
#  
#  Increasing the scale makes recognition-lengths longer, a scale
#  of 100 will force full-length recognition.
#  
#  Set the skew (-k) or N-pct (-x) to 100 or 0 to turn it off
