#!/bin/bash

out=$1
shift

if [ $# -gt 1 ] ; then
    samtools merge $out $@
elif [ $# -eq 1 ] ; then
    cp $1 $out
else
    echo 'Need to specify some arguments, you know' >&2
    echo 'Usage: $0 out.bam in.bam [in2.bam ...]' >&2
    exit 1
fi
