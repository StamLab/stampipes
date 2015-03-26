SAMPLE_SIZE=$1
FAI=$2
INBAM=$3
OUTBAM=$4

N_READS=`samtools view -c ${INBAM}`

if [ ${SAMPLE_SIZE} -gt ${N_READS} ]
then
    if [ ! -e $OUTBAM ]
    then
        ln -s $INBAM $OUTBAM
    fi
else
    samtools view ${INBAM} | 
		random-lines -n${SAMPLE_SIZE} -N${N_READS} | \
		samtools view -bS -t ${FAI} - > \
		${OUTBAM}
fi
