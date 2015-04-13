SAMPLE_NAME=$1
INBAM=$2
OUTPUT=$3

# Delete output if it already exists; do not want to append to existing output
if [ -e $OUTPUT ]; then
    rm $OUTPUT
fi

echo "Calculating total/pf/qc counts"
if [ -e ${SAMPLE_NAME}_R1_001.fastq.gz ]; then
zcat ${SAMPLE_NAME}_R?_???.fastq.gz \
     | awk 'BEGIN{
     filter = 0;
   }{ \
   if ( FNR % 4 == 1 && substr($2, 3, 1) == "Y" ) \
     { filter+=1 }
   }\
   END { \
     TOTAL = NR / 4; \
     print "total", TOTAL ; \
     print "pf", TOTAL - filter ; \
     print "qc", filter ; \
   }' >> $OUTPUT
fi

echo "Calculate tags trimmed"
# count is in read pairs, double for accuracy in LIMS display % wise
TRIMMED=`find . -maxdepth 1 -name "$SAMPLE_NAME*trimstats.txt" | xargs awk 'BEGIN { total=0 } { match($0, /Total read-pairs trimmed: ([0-9]+)/, a); total=total+a[1];} END{ print total * 2 }'`
echo -e "adapter-trimmed\t$TRIMMED" >> $OUTPUT

echo "Creating bam counts"
python $STAMPIPES/scripts/bwa/bamcounts.py $INBAM $OUTPUT
