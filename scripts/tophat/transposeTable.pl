#!/usr/bin/perl -w
use strict;

# Script to transpose a tab-delimited text file
my @arrayRefs = ();
my $maxCol = 0;
while (<>) {
    chomp;
    my @array = split /\t/;
    push @arrayRefs, \@array;
    my $cols = @array;
    if ($cols > $maxCol) {
        $maxCol = $cols;
    }
}

my $lines = @arrayRefs;
for (my $i = 0; $i < $maxCol; ++$i ) {
    my $line = "";
    for (my $j = 0; $j < $lines; ++$j ) {
        my $arrayRef = $arrayRefs[$j];
        my $value = $arrayRef->[$i];
        unless (defined($value)) {
            $value = "";
        }
        if ($j > 0) {
            $line .= "\t";
        }
        $line .= $value;
    }
    print "$line\n";
}
