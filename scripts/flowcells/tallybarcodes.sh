#!/bin/bash

# Usage: $0 1.fastq.gz [ 2.fastq.gz ... ]
# Barcodes are the are the last colon-separated field of every 4th line
# sort | uniq -c | sort is slow, so we use a perl implementation
zcat "$@" \
  | sed -n '1~4p' \
  | awk -F: '{print $NF}' \
  | perl -ne '$seen{$_}++ ; END{ for $k (sort {$seen{$b} <=> $seen{$a} } keys %seen){ print "$seen{$k}\t$k"; } };'
