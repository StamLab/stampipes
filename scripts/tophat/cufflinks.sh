#!/bin/bash

# This is to be run within a cluster script
# Quantitates against the GTF, does not assemble transcripts de-novo,
# that would be a special case for deep paired-end datasets
# Should I add an option here?
echo "$0 $1 $2 $3 $4 $5 on "`hostname`" at "`date`

if [ "$#" -lt "7" ] ; then
    echo "Usage:  $0   alignment.bam genome library-type reference_annotations.gtf output.dir number.of.threads mode=[GTF|GTF-guide]"
    exit 100;
fi
bamfile=$1
refseq=$2
libraryType=$3
gtf=$4
outputdir=$5
numthreads=$6
gtfmode=$7

if ! [ -d "$outputdir" ] ; then
    mkdir -p $outputdir
fi

refseqfasta=$REF_DIR/$refseq.fa

#if [ "$TMPDIR_t" == "" ] ; then
#    TMPDIR_t=/tmp
#fi
#echo "TMPDIR_t=$TMPDIR_t"

if ! [ -d "$outputdir" ] ; then
    mkdir -p $outputdir
fi

hostname
date


optparams=" --multi-read-correct        "
if [ "$refseq" != "hg19" ] ; then
    # This crashes on Gencode annotations only
    optparams="$optparams -b $refseqfasta "
fi
        #-G $gtf \
cufflinks -v --no-update-check -p $numthreads -o $outputdir \
        --$gtfmode $gtf \
        --library-type $libraryType \
        --multi-read-correct        \
        $bamfile

if [ "$?" -ne "0" ]; then
    exit 100
else
    exit 0
fi

#  cufflinks v1.2.1
#  linked against Boost version 104000
#  -----------------------------
#  Usage:   cufflinks [options] <hits.sam>
#  General Options:
#    -o/--output-dir              write all output files to this directory              [ default:     ./ ]
#    -p/--num-threads             number of threads used during analysis                [ default:      1 ]
#    --seed                       value of random number generator seed                 [ default:      0 ]
#    -G/--GTF                     quantitate against reference transcript annotations                      
#    -g/--GTF-guide               use reference transcript annotation to guide assembly                   
#    -M/--mask-file               ignore all alignment within transcripts in this file                     
#    -b/--frag-bias-correct       use bias correction - reference fasta required        [ default:   NULL ]
#    -u/--multi-read-correct      use 'rescue method' for multi-reads (more accurate)   [ default:  FALSE ]
#    --library-type               library prep used for input reads                     [ default:  below ]
#  
#  Advanced Abundance Estimation Options:
#    -m/--frag-len-mean           average fragment length (unpaired reads only)         [ default:    200 ]
#    -s/--frag-len-std-dev        fragment length std deviation (unpaired reads only)   [ default:     80 ]
#    --upper-quartile-norm        use upper-quartile normalization                      [ default:  FALSE ]
#    --max-mle-iterations         maximum iterations allowed for MLE calculation        [ default:   5000 ]
#    --num-importance-samples     number of importance samples for MAP restimation      [ default:   1000 ]
#    --compatible-hits-norm       count hits compatible with reference RNAs only        [ default:  FALSE ]
#    --total-hits-norm            count all hits for normalization                      [ default:  TRUE  ]
#  
#  Advanced Assembly Options:
#    -L/--label                   assembled transcripts have this ID prefix             [ default:   CUFF ]
#    -F/--min-isoform-fraction    suppress transcripts below this abundance level       [ default:   0.10 ]
#    -j/--pre-mrna-fraction       suppress intra-intronic transcripts below this level  [ default:   0.15 ]
#    -I/--max-intron-length       ignore alignments with gaps longer than this          [ default: 300000 ]
#    -a/--junc-alpha              alpha for junction binomial test filter               [ default:  0.001 ]
#    -A/--small-anchor-fraction   percent read overhang taken as 'suspiciously small'   [ default:   0.09 ]
#    --min-frags-per-transfrag    minimum number of fragments needed for new transfrags [ default:     10 ]
#    --overhang-tolerance         number of terminal exon bp to tolerate in introns     [ default:      8 ]
#    --max-bundle-length          maximum genomic length allowed for a given bundle     [ default:3500000 ]
#    --max-bundle-frags           maximum fragments allowed in a bundle before skipping [ default: 500000 ]
#    --min-intron-length          minimum intron size allowed in genome                 [ default:     50 ]
#    --trim-3-avgcov-thresh       minimum avg coverage required to attempt 3' trimming  [ default:     10 ]
#    --trim-3-dropoff-frac        fraction of avg coverage below which to trim 3' end   [ default:    0.1 ]
#  
#  Advanced Reference Annotation Guided Assembly Options:
#    --no-faux-reads              disable tiling by faux reads                          [ default:  FALSE ]
#    --3-overhang-tolerance       overhang allowed on 3' end when merging with reference[ default:    600 ]
#    --intron-overhang-tolerance  overhang allowed inside reference intron when merging [ default:     30 ]
#  
#  Advanced Program Behavior Options:
#    -v/--verbose                 log-friendly verbose processing (no progress bar)     [ default:  FALSE ]
#    -q/--quiet                   log-friendly quiet processing (no progress bar)       [ default:  FALSE ]
#    --no-update-check            do not contact server to check for update availability[ default:  FALSE ]
#  
#  Supported library types:
#  	ff-firststrand
#  	ff-secondstrand
#  	ff-unstranded
#  	fr-firststrand
#  	fr-secondstrand
#  	fr-unstranded (default)
#  	transfrags

