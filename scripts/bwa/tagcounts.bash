SAMPLE_NAME=$1
INBAM=$2
OUTPUT=$3

# Delete output if it already exists; do not want to append to existing output
if [ -e $OUTPUT ]; then
    rm $OUTPUT
fi

python $STAMPIPES/scripts/bwa/bamcounts.py $INBAM $OUTPUT

# Add in total/PF/QC counts
# TODO: Can probably achieve TOTAL/PF/QC in a single pass of awk
echo "Calculating total count"
TOTAL=`zcat ${SAMPLE_NAME}_R?_???.fastq.gz | wc -l | awk '{print $1 / 4; }'`
echo -e "total\t$TOTAL" >> $OUTPUT

echo "Calculating pf count"
PF=`zcat ${SAMPLE_NAME}_R?_???.fastq.gz | \
  awk '{if (substr($2, 3, 1) == "N") {f=0;print $1} else if (substr($2, 3, 1) == "Y") {f=1} else if ( f == 0) {print $1 } }' | \
  wc -l | awk '{print $1 / 4; }'`
echo -e "pf\t$PF" >> $OUTPUT

echo "Calculating qc count"
QC=`awk -v TOTAL=$TOTAL -v PF=$PF 'BEGIN{ print TOTAL - PF; }'`
echo -e "qc\t$QC" >> $OUTPUT
