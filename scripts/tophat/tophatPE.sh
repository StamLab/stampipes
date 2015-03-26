#!/bin/bash

# This is to be run within a cluster script
# Aligns a single pair of reads files

if [ "$#" -lt "6" ] ; then
    echo "Usage:  $0   readsF.fastq.gz readsR.fastq.gz genome library-type annot.gtf output.dir number.of.threads" 1>&2
    exit 100
fi
readsF=$1
readsR=$2
refseq=$3
libraryType=$4
gtf=$5
outputdir=$6
numthreads=$7

if ! [ -d "$outputdir" ] ; then
    mkdir -p $outputdir
fi

bowtie_index=$REF_DIR/$refseq

if ! [ -s "$gtf" ] ; then
    echo "Failed to find $gtf"
    exit 1
fi


# No idea at all what this should be
mate_inner_dist=200 
# Need to set this larger than default
mate_std_dev=300

hostname
date

opt_known_junctions=" --GTF $gtf "
if [ "$IS_MICROBE_RNA" == "1" ] ; then
    echo "Omitting known junctions for bacterial RNAseq, might as well be using bowtie or anything else" 1>&2
    opt_known_junctions=" "
fi

opt_special=""
echo -e "USE_TMERCER_PARAMS=($USE_TMERCER_PARAMS)"
if [ "$USE_TMERCER_PARAMS" == "1" ] ; then
    echo "Using Tim Mercer's mapping parameters for TopHat"
    # reducing -g/--max-multihits from default 20 to 10 
    # reducing --segment-length from default 25 to 18
    # reducing --segment-mismatches from default 2 to 0
    opt_special=" -g 10 --segment-length 18 --segment-mismatches 0 "
fi

#options=" --GTF $gtf -r $mate_inner_dist --mate-std-dev $mate_std_dev --tmp-dir $TMPDIR_t --library-type $libraryType "
options=" $opt_known_junctions -r $mate_inner_dist --mate-std-dev $mate_std_dev --library-type $libraryType --no-coverage-search"
reads=" $readsF $readsR "
echo -e "tophat $options $opt_special -p $numthreads -o $outputdir $bowtie_index $reads"
         tophat $options $opt_special -p $numthreads -o $outputdir $bowtie_index $reads
exitstatus=$?

if [ "$exitstatus" -ne "0" ]; then
    echo "$0:  TopHat failed with exit code ($exitstatus)" 1>&2
    exit 100
else
    resultfile="$outputdir/accepted_hits.bam"
    if ! [ -s "$resultfile" ] ; then
        echo "$0:  Failed to create $resultfile" 1>&2
        exit 100
    else
        exit 0
    fi
fi

#  tophat: 
#  TopHat maps short sequences from spliced transcripts to whole genomes.
#  
#  Usage:
#      tophat [options] <bowtie_index> <reads1[,reads2,...]> [reads1[,reads2,...]] \
#                                      [quals1,[quals2,...]] [quals1[,quals2,...]]
#  
#  Options:
#      -v/--version
#      -o/--output-dir                <string>    [ default: ./tophat_out     ]
#      -a/--min-anchor                <int>       [ default: 8                ]
#      -m/--splice-mismatches         <0-2>       [ default: 0                ]
#      -i/--min-intron-length         <int>       [ default: 50               ]
#      -I/--max-intron-length         <int>       [ default: 500000           ]
#      -g/--max-multihits             <int>       [ default: 20               ]
#      -F/--min-isoform-fraction      <float>     [ default: 0.15             ]
#      --max-insertion-length         <int>       [ default: 3                ]
#      --max-deletion-length          <int>       [ default: 3                ]
#      --solexa-quals
#      --solexa1.3-quals                          (same as phred64-quals)
#      --phred64-quals                            (same as solexa1.3-quals)
#      -Q/--quals
#      --integer-quals
#      -C/--color                                 (Solid - color space)
#      --color-out
#      --library-type                 <string>    (fr-unstranded, fr-firststrand,
#                                                  fr-secondstrand)
#      -p/--num-threads               <int>       [ default: 1                ]
#      -G/--GTF                       <filename>
#      -j/--raw-juncs                 <filename>
#      --insertions                   <filename>
#      --deletions                    <filename>
#      -r/--mate-inner-dist           <int>
#      --mate-std-dev                 <int>       [ default: 20               ]
#      --no-novel-juncs
#      --no-novel-indels
#      --no-gtf-juncs
#      --no-coverage-search
#      --coverage-search
#      --no-closure-search
#      --closure-search
#      --microexon-search
#      --butterfly-search
#      --no-butterfly-search
#      --keep-tmp
#      --tmp-dir                      <dirname>   [ default: <output_dir>/tmp ]
#      -z/--zpacker                   <program>   [ default: gzip             ]
#      --no-unmapped-fifo
#  
#  Advanced Options:
#      --segment-mismatches           <int>       [ default: 2                ]
#      --segment-length               <int>       [ default: 25               ]
#      --bowtie-n                                 [ default: bowtie -v        ]
#      --min-closure-exon             <int>       [ default: 100              ]
#      --min-closure-intron           <int>       [ default: 50               ]
#      --max-closure-intron           <int>       [ default: 5000             ]
#      --min-coverage-intron          <int>       [ default: 50               ]
#      --max-coverage-intron          <int>       [ default: 20000            ]
#      --min-segment-intron           <int>       [ default: 50               ]
#      --max-segment-intron           <int>       [ default: 500000           ]
#      --no-sort-bam                              [Output BAM is not coordinate-sorted]
#      --no-convert-bam                           [Do not convert to bam format.
#                                                  Output is <output_dir>accepted_hit.sam.
#                                                  Implies --no-sort-bam.]
#  
#  SAM Header Options (for embedding sequencing run metadata in output):
#      --rg-id                        <string>    (read group ID)
#      --rg-sample                    <string>    (sample ID)
#      --rg-library                   <string>    (library ID)
#      --rg-description               <string>    (descriptive string, no tabs allowed)
#      --rg-platform-unit             <string>    (e.g Illumina lane ID)
#      --rg-center                    <string>    (sequencing center name)
#      --rg-date                      <string>    (ISO 8601 date of the sequencing run)
#      --rg-platform                  <string>    (Sequencing platform descriptor)
#  
#      for detailed help see http://tophat.cbcb.umd.edu/manual.html




