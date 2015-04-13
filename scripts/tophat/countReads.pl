#!/usr/bin/perl -w
use strict;

# add up lines in fastq files, print name key value

if (@ARGV < 3) {
    die "Usage:  $0 name key [reads directory]\n";
}

my ($name,$key,$readsDir) = @ARGV;

unless (-d $readsDir) {
    die "($readsDir) directory not found\n";
}

my $sum = 0;
foreach my $file (glob("$readsDir/${name}_R?_???.fastq.gz")) {
    chomp( my $lines = `zcat $file | wc -l` );
    $sum += ($lines / 4);
}

print "$name\t$key\t$sum\n";

