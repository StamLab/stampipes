#!/usr/bin/env bash
# make_rRNA.sh
# Kamil Slowikowski
# December 12, 2014
#
# 1. Download chromosome sizes from UCSC if needed.
# 2. Make an interval_list file suitable for CollectRnaSeqMetrics.jar.
#
# Gencode v19 genes:
#
#   ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
#
# Picard Tools CollectRnaSeqMetrics.jar:
#
#   https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics

# 1. Chromosome sizes from the UCSC genome browser ---------------------------

chrom_sizes=hg19.chrom.sizes

if [[ ! -s $chrom_sizes ]]
then
    mysql -N --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
        "SELECT chrom,size FROM chromInfo ORDER BY size DESC;" hg19 \
    > $chrom_sizes
fi

# 2. rRNA interval_list file -------------------------------------------------

# Genes from Gencode.
genes=gencode.v19.annotation.gtf

# Output file suitable for Picard CollectRnaSeqMetrics.jar.
rRNA=gencode.v19.rRNA.interval_list

# Sequence names and lengths. (Must be tab-delimited.)
perl -lane 'print "\@SQ\tSN:$F[0]\tLN:$F[1]\tAS:hg19"' $chrom_sizes | \
    grep -v _ \
>> $rRNA

# Intervals for rRNA transcripts.
grep 'gene_type "rRNA"' $genes | \
    awk '$3 == "transcript"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '
        /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
        print join "\t", (@F[0,1,2,3], $1)
    ' | \
    sort -k1V -k2n -k3n \
>> $rRNA