SAMPLE_SIZE=$1
FAI=$2
INFILE=$3
OUTFILE=$4

N_READS=`wc -l ${INFILE} | cut -f1 -d" "`

if [ ${SAMPLE_SIZE} -gt ${N_READS} ]
then
    awk -f $STAMPIPES/awk/splitreadpairs.awk $INFILE | \
    samtools view -bSt $FAI - > $OUTFILE
else
    cat ${INFILE} |
        random-lines -n${SAMPLE_SIZE} -N${N_READS} | \
        awk -f $STAMPIPES/awk/splitreadpairs.awk | \
        samtools view -bSt $FAI - > \
        ${OUTFILE}
fi
