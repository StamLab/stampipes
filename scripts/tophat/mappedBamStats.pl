#!/usr/bin/perl -w
use strict;

if (@ARGV < 1) {
    die "Usage:  $0 file.bam [name to override first dot-separated token in bam file basename]\n";
}
my ($bamfile,$name) = @ARGV;
unless (defined($name)) {
    chomp( $name=`basename $bamfile | cut -f1 -d .` );
}

my $mappedTotal = 0;
my $mappedChrM = 0;
my $duplicates = 0;
#my $readcmd = "samtools view -X $bamfile";
my $readcmd = "samtools view $bamfile";
open IN, "$readcmd |" or die "Failed to $readcmd\n";
while (<IN>) {
    chomp;
    my ($readname,$flags,$chrom,$position,$mappingQuality,@other) = split /\t/;
    my $isSecondaryAlignment = ($flags & 0x0100);
    if ((! $isSecondaryAlignment) and ($chrom ne "*")) {
        ++$mappedTotal;
        if ("chrM" eq $chrom) {
            ++$mappedChrM;
        }
        if ($flags & 0x0400) {
            ++$duplicates;
        }
    }
    # Okay I've figured this out.  Ignore TopHat mapping quality.

    # notes on TopHat mapping quality, unique alignments
    # TopHat alignment output BAM files contain multiple entries for some reads. Sometimes one position is better, sometimes it is completely ambiguous.
    # Totals from five NuGEN libraries:
    # 
    # primary alignments (exactly one per input read that maps at least one place):
    # mapQ       count
    #   0       478,216
    #   1     1,183,594
    #   3     1,485,285   * same as secondary with score 3
    # 255    13,981,891
    # 
    # secondary alignments:
    #   0 3,504,028
    #   1 2,629,314
    #   3 1,485,285   *
    # So mapQ 3 must mean exactly two equivalent positions? 
    # And mapQ 1 and 0 must represent larger numbers of possible positions. 
    # Information from seqanswers forum:
    # mapQ 255 = only one position found. TopHat found no particular worse position for comparison, so sets quality to undefined. 
    # Downstream, Cufflinks apparently has a smart way to deal with multi-mapping reads, using them for transcript assembly (avoiding gaps between unique coverage) but not for RPKM. 
    # 
    # Bottom line: I plan to include multi-mapping reads in coverage plots, and use Cufflinks to calculate RPKMs.



    # warn "2ndary($isSecondaryAlignment) mapQ($mappingQuality)\n";
    #if ($isSecondaryAlignment xor (0 == $mappingQuality)) {
    #    warn "inconsistent 2ndary($isSecondaryAlignment) mapQ($mappingQuality)\n";
    #}
    #foreach my $flagchar (split //, $flags) { ++$flagcounts{$flagchar}; }
}
close IN;

print "$name\tmapped_reads\t$mappedTotal\n";
print "$name\tmapped_chrM\t$mappedChrM\n";
print "$name\tduplicates\t$duplicates\n";

my $flagsdoc=<<EOFl;
Flag    Chr     Description
0x0001  p       the read is paired in sequencing
0x0002  P       the read is mapped in a proper pair
0x0004  u       the query sequence itself is unmapped
0x0008  U       the mate is unmapped
0x0010  r       strand of the query (1 for reverse)
0x0020  R       strand of the mate
0x0040  1       the read is the first read in a pair
0x0080  2       the read is the second read in a pair
0x0100  s       the alignment is not primary
0x0200  f       the read fails platform/vendor quality checks
0x0400  d       the read is either a PCR or an optical duplicate
EOFl

#   my %flagcounts = ();
#   foreach my $flagchar qw( p P u U r R 1 2 s f d ) {
#       $flagcounts{$flagchar} = 0;
#   }
#   my %flagdesc = (
#   p =>  "the read is paired in sequencing",
#   P =>  "the read is mapped in a proper pair",
#   u =>  "the query sequence itself is unmapped",
#   U =>  "the mate is unmapped",
#   r =>  "strand of the query (1 for reverse)",
#   R =>   "strand of the mate",
#   1 =>   "the read is the first read in a pair",
#   2 =>   "the read is the second read in a pair",
#   s =>   "the alignment is not primary",
#   f =>   "the read fails platform/vendor quality checks",
#   d =>   "the read is either a PCR or an optical duplicate"
#   );
#   
#   foreach my $flagchar (sort keys %flagcounts) {
#       my $count = $flagcounts{$flagchar};
#       if ($count > 0) {
#           my $desc = $flagdesc{$flagchar};
#           $desc = join("_", (split / /, $desc) );
#           print "$name\tflag-$flagchar-$desc\t$count\n";
#       }
#   }
