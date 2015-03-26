#!/usr/bin/perl -w
use strict;

my $prevChrom = undef;
my $prevMin0 = undef;
my $prevMax1 = undef;
my $prevScore = undef;
while (<STDIN>) {
    chomp;
    my ($chrom,$min0,$max1,$score) = split;
    if (defined( $prevChrom ) and ($prevChrom eq $chrom) and ($prevMax1 == $min0) and ($prevScore == $score)) {
        $prevMax1 = $max1; # combine these two adjacent intervals with the same score
    } else {
        printLast();
        ($prevChrom,$prevMin0,$prevMax1,$prevScore) = ($chrom,$min0,$max1,$score);
    }
}
printLast();

# print and clear last data 
sub printLast {
    if (defined($prevChrom)) {
        print "$prevChrom\t$prevMin0\t$prevMax1\t$prevScore\n";
    }
}
