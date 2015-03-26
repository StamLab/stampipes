#!/usr/bin/perl -w
use strict;

# Just add placeholder commas to any space-delimited field thats all digits
# Also turns any whitespace into tabs, like it or not.

while (<>) {
    chomp;
    my $result = "";
    foreach my $field (split) {
        if ($field =~ /^([\d]+)$/) {
            $field = encomma($field);
        }
        if (length($result) > 0) {
            $result .= "\t";
        }
        $result .= $field;
    }
    $result .= "\n";
    print $result;
}

sub encomma {
    my $orig = shift;
    my $remain = "$orig";
    my $result = "";
    while (length($remain) > 3) {
        my $chunk = substr( $remain, length($remain) - 3 );
        $remain = substr( $remain, 0, length($remain) - 3 );
        if (length($result) > 0) {
            $result = ",$result";
        }
        $result = "$chunk$result";
        #warn "$remain | $result\n";
    }
    if ((length($result) > 0) and (length($remain)>0)) {
        $result = ",$result";
    }
    $result = "$remain$result";
    #warn "encomma($orig)=($result)\n";
    return $result;
}

