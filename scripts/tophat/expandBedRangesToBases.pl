#!/usr/bin/perl -w
use strict;

#my $regions = 0;
while (<>) {
    chomp;
    my ($chrom,$min,$max,@whatever) = split /\t/;
    #my $tmpname = "region$regions";
    #++$regions;
    for (my $offset = 0; $min + $offset < $max; ++$offset ) {
        my $min0 = $min + $offset;
        my $max1 = $min0 + 1;
        #print "$chrom\t$min0\t$max1\t$tmpname\t$offset\n";
        print "$chrom\t$min0\t$max1\n"; # bed3
    }
}

