#!/usr/bin/perl -w
use strict;

#die "ERROR:  $0 replaced by splitCoverageByTranscriptStrand.pl\n";

if (@ARGV < 3) {
    die "Usage:  $0 input.bam output_prefix genome_or_other_suffix [optional chromosome]\n" .
        "       writes 3 files: \$output_prefix.[pos|neg|all].\$genome.starch\n";
}
my ($inputBam, $outputPrefix, $outputLabel2, $chromosome) = @ARGV;

my $bam2bed = $ENV{'SCRIPT_DIR'} . "/bam2bed";

open OUTALL, "| starch - > $outputPrefix.all.$outputLabel2.starch" or die "$!";
open OUTPOS, "| starch - > $outputPrefix.pos.$outputLabel2.starch" or die "$!";
open OUTNEG, "| starch - > $outputPrefix.neg.$outputLabel2.starch" or die "$!";

my $cmd = "$bam2bed --split < $inputBam "; 
if (defined($chromosome)) {
    $cmd = "samtools view -u $inputBam $chromosome | $bam2bed --split "; 
}
# This could maybe be changed to samtools view?
open IN, "$cmd |" or die "$!";
while (<IN>) {
    chomp;
    my ($chrom,$min0,$max1,$qname,$flags,$strand,@ignored) = split /\t/;
    #    - FLAG                      <-->   score (5th column)
    #    - 16 & FLAG                 <-->   strand (6th column)
    my $bed3 = "$chrom\t$min0\t$max1\n";
    print OUTALL $bed3;
    my $reverse = $flags & 0x10;
    my $isR1 = $flags & 0x40;
    my $isR2 = $flags & 0x80;
    if ($reverse) {
        # "-" strand, opposite of reference
        if ($isR1) {
            print OUTNEG $bed3;
        } elsif ($isR2) {
            print OUTPOS $bed3; # double negative
        } else {
            die $_;
        }
    } else {
        # "+" strand, same as reference
        if ($isR1) {
            print OUTPOS $bed3;
        } elsif ($isR2) {
            print OUTNEG $bed3;
        } else {
            die $_;
        }
    }
    #warn "SANITY?\t$strand\t$reverse\t$isR1\t$isR2\n";
}
close IN;

close OUTALL;
close OUTPOS;
close OUTNEG;



#
#   if [ "$#" -lt "2" ] ; then
#       echo "Usage:  $0 file.bam genome_build"
#       exit 1;
#   fi
#
#   bam=$1
#   # refseq chrom.sizes required by bedGraphToBigWig
#   refseq=$2
#
#   name=`basename $bam | cut -f1 -d .`
#   destdir=`dirname $bam`
#   gcovpos=$destdir/coverage.$name.pos.$refseq.bed
#   gcovneg=$destdir/coverage.$name.neg.$refseq.bed
#   gcovall=$destdir/coverage.$name.all.$refseq.bed
#
#   # First see if this has already been done
#   #for file in $gcovpos $gcovneg ; do
#   for file in $gcovall ; do
#       subname=${file%%.bed*}
#       bigwig=$subname.bigWig
#       if [ -s "$bigwig" ] ; then
#           echo "Already have $bigwig" 1>&2
#           exit
#       fi
#   done
#
#   # Create the coverage bed file for each strand
#   if ! [ -f "$gcovpos" ] ; then
#       destpos=$destdir/readspans.$name.pos.$refseq.bed3
#       destneg=$destdir/readspans.$name.neg.$refseq.bed3
#       destall=$destdir/readspans.$name.all.$refseq.bed3
#       if ! [ -f "$destpos" ] ; then
#           echo -e "$name $refseq"
#           bamToBed -split -i $bam \
#               | $SCRIPT_DIR/bbms.pl \
#               | awk "BEGIN{OFS=\"\\t\"}{ print \$1,\$2,\$3 > \"$destall\"; if (\$6==\"+\") { print \$1,\$2,\$3 > \"$destpos\" } else if (\$6==\"-\") { print \$1,\$2,\$3 > \"$destneg\" }; }"
#       fi
#       $SCRIPT_DIR/singleBedFileBaseCoverage.sh $destpos | compressBed4.pl > $gcovpos
#       #cat $gcovpos | awk '{print $3-$2}' | add2
#       singleBedFileBaseCoverage.sh $destneg | compressBed4.pl > $gcovneg
#       #cat $gcovneg | awk '{print $3-$2}' | add2
#       singleBedFileBaseCoverage.sh $destall | compressBed4.pl > $gcovall
#       #cat $gcovall | awk '{print $3-$2}' | add2
#       if [ -s "$gcovpos" ] ; then
#           # cleanup
#           rm $destpos $destneg $destall
#       fi
#   fi
#
#   chromSizes=$SCRIPT_DIR/refseq/$refseq/$refseq.chrom.sizes
#   if ! [ -f "$chromSizes" ] ; then
#       chromSizes=$SCRIPT_DIR/refseq/$refseq/chrom_sizes.txt
#       if ! [ -f "$chromSizes" ] ; then
#           echo "Failed to find $chromSizes " 1>&2
#           exit 1
#       fi
#   fi
#
#   # Convert to bigWig
#   for file in $gcovpos $gcovneg $gcovall ; do
#       name=${file%%.bed*}
#       bigwig=$name.bigWig
#       if ! [ -f "$bigwig" ] ; then
#           $SCRIPT_DIR/bedGraphToBigWig $file $chromSizes $bigwig
#       fi
#       if [ -s "$bigwig" ] ; then
#           # cleanup
#           rm $file
#       fi
#   done
#
#
#  Usage:
#    bam2bed [ --help ] [ --keep-header ] [ --split ] [ --all-reads ] [ --do-not-sort | --max-mem <value> (--sort-tmpdir <dir>) ] < foo.bam
#
#  Options:
#    --help                 Print this help message and exit
#    --keep-header          Preserve header section as pseudo-BED elements
#    --split                Split reads with 'N' CIGAR operations into separate BED elements
#    --all-reads            Include both unmapped and mapped reads in output
#    --do-not-sort          Do not sort converted data with BEDOPS sort-bed
#    --custom-tags <value>  Add a comma-separated list of custom SAM tags
#    --max-mem <value>      Sets aside <value> memory for sorting BED output. For example,
#                           <value> can be 8G, 8000M or 8000000000 to specify 8 GB of memory
#                           (default: 2G)
#    --sort-tmpdir <dir>    Optionally sets <dir> as temporary directory for sort data, when
#                           used in conjunction with --max-mem <value>. For example, <dir> can
#                           be $PWD to store intermediate sort data in the current working
#                           directory, in place of the host operating system default
#                           temporary directory.
#
#  About:
#    This script converts 0-based, half-open [a-1, b) binary BAM data from standard
#    input into 0-based, half-open [a-1, b) extended BED, which is sorted and sent
#    to standard output.
#
#    We process BAM columns from mappable reads (as described by specifications at
#    http://samtools.sourceforge.net/SAM1.pdf) converting them into the first six
#    UCSC BED columns as follows:
#
#    - RNAME                     <-->   chromosome (1st column)
#    - POS - 1                   <-->   start (2nd column)
#    - POS + length(CIGAR) - 1   <-->   stop (3rd column)
#    - QNAME                     <-->   id (4th column)
#    - FLAG                      <-->   score (5th column)
#    - 16 & FLAG                 <-->   strand (6th column)
#
#    The remaining BAM columns are mapped intact, in order, to adjacent BED columns:
#
#    - MAPQ
#    - CIGAR
#    - RNEXT
#    - PNEXT
#    - TLEN
#    - SEQ
#    - QUAL
#
#    Use the --keep-header option if you would like to preserve the SAM header section
#    as pseudo-BED elements; otherwise, these are stripped from output.
#
#    Because we have mapped all columns, we can translate converted BED data back to
#    headerless SAM reads with a simple awk statement or other script that calculates
#    1-based coordinates and permutes columns.
#
#    By default, we only process mapped reads. If you also want to convert unmapped
#    reads, add the --all-reads option.
#
#    In the case of RNA-seq data with skipped regions ('N' components in the read's
#    CIGAR string), the --split option will split the read into two or more separate
#    BED elements.
#
#    This script also validates the CIGAR strings in a sequencing dataset, in the
#    course of converting to BED.
#
#    Example usage:
#
#    $ bam2bed < foo.bam > sorted-foo.bam.bed
#
#    Because we make no assumptions about the sort order of your input, we apply the
#    sort-bed application to output to generate lexicographically-sorted BED data.
#
#    If you want to skip sorting, use the --do-not-sort option:
#
#    $ bam2bed --do-not-sort < foo.bam > unsorted-foo.bam.bed
#
#    This option is *not* recommended, however, as other BEDOPS tools require sorted
#    inputs to process data correctly.
