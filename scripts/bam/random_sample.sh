#!/usr/bin/env bash
# 
USAGE_MSG="Usage: $0 in.bam out.bam numPairs [pairTotal]"


die(){
  echo $@ >&2
  exit 1
}

if [ $# -lt 3 ] ; then
  die "$USAGE"
fi

inbam="$1"
outbam="$2"
numPairs="$3"
pairTotal="$4"

# If total size of file is not specified, read from samtools
if [ -z "$pairTotal" ] ; then
  if [ -s "$inbam.bai" ] ;then
    # If index exists, use idxstats
    readTotal=$(samtools idxstats $inbam | awk '{sum+=$3}END{print sum}')
    pairTotal=$((readTotal / 2))
  else
    # Otherwise use view to count them
    pairTotal=$(samtools view -c -f 64 $inbam)
  fi

fi

echo "Random sampling started at: `date +%D%t%T.%N`"
echo "Sampling $numPairs read-pairs from $inbam to $outbam (from $pairTotal total)..."

python3 $STAMPIPES/scripts/bam/random_reads.py "$inbam" "$outbam" "$pairTotal" "$numPairs"

echo "Random sampling finished at: `date +%D%t%T.%N`"
