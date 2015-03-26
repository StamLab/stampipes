#!/usr/bin/perl -w
use strict;

# add up lines in files, print name key value

if (@ARGV < 3) {
    die "Usage:  $0 name key file1 [file2]...\n";
}

my ($name,$key,@files) = @ARGV;

my $sum = 0;
foreach my $file (@files) {
    unless (-f $file) {
        die "Failed to find $file\n";
    }
    chomp( my $lines = `zcat -f $file | wc -l` );
    $sum += $lines;
}

print "$name\t$key\t$sum\n";
