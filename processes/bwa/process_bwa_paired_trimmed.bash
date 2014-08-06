FINAL_BAM=${SAMPLE_NAME}.sorted.bam
UNIQUES_BAM=${SAMPLE_NAME}.uniques.sorted.bam

if [ ! -e ${SAMPLE_NAME}_R1_fastqc -o ! -e ${SAMPLE_NAME}_R2_fastqc ]; then
qsub -N "fq${SAMPLE_NAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`
  
  make -f $STAMPIPES/makefiles/fastqc.mk
__SCRIPT__
fi

NUMBER_FASTQ_FILES=`find . -name "${SAMPLE_NAME}_R1_???.fastq.gz" | wc -l`
FASTQ_PAIR_HOLDS=""
FASTQ_PAIR_BAMS=""

for filenum in $(seq -f "%03g" 1 $NUMBER_FASTQ_FILES)
do
  NAME=".aln${SAMPLE_NAME}_${filenum}_${FLOWCELL}"
  BAMFILE="${SAMPLE_NAME}_${filenum}.sorted.bam"
  
if [ ! -e $BAMFILE -a ! -e ${FINAL_BAM} ]; then
    
qsub -l h_data=5650M -N ${NAME} -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`
  
  make -f $STAMPIPES/makefiles/bwa/bwa_paired_trimmed_make.mk \
    FASTQ1_FILE=${SAMPLE_NAME}_R1_${filenum}.fastq.gz \
    FASTQ2_FILE=${SAMPLE_NAME}_R2_${filenum}.fastq.gz \
    OUTBAM=${BAMFILE} \
    TRIMSTATS=${SAMPLE_NAME}_${filenum}.trimstats.txt \
    ADAPTERFILE=$STAMPIPES/data/adapters/default.adapters
__SCRIPT__

  # Only hold on alignments that are being run
  FASTQ_PAIR_HOLDS="$FASTQ_PAIR_HOLDS,$NAME"
fi

# Need to keep track of these even if they have already finshed
# for proper merging
FASTQ_PAIR_BAMS="${BAMFILE} ${FASTQ_PAIR_BAMS}"

done

if [ -n "$FASTQ_PAIR_HOLDS" ]; then
    HOLD="-hold_jid $FASTQ_PAIR_HOLDS"
fi

if [ ! -e ${FINAL_BAM} -o ! -e ${UNIQUES_BAM} ]; then

PROCESS_HOLD="-hold_jid .pb${SAMPLE_NAME}_${FLOWCELL}"
    
qsub ${HOLD} -N ".pb${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`
  
  if [ ! -e ${FINAL_BAM} ]; then
  if [ "$NUMBER_FASTQ_FILES" -eq "1" ]
  then
    mv ${SAMPLE_NAME}_001.sorted.bam ${SAMPLE_NAME}.sorted.bam
  else
    samtools merge ${FINAL_BAM} ${FASTQ_PAIR_BAMS}
  fi
  fi
  
  make -f $STAMPIPES/makefiles/bwa/process_paired_bam.mk

  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
    -t ${LIMS_API_TOKEN} \
    -f ${FLOWCELL} \
    --alignment_id ${ALIGNMENT_ID} \
    --flowcell_lane_id ${FLOWCELL_LANE_ID} \
    --insertsfile ${SAMPLE_NAME}.CollectInsertSizeMetrics.picard
  
  if [ "$NUMBER_FASTQ_FILES" -gt "1" ]
  then
    rm $FASTQ_PAIR_BAMS
  fi
__SCRIPT__

fi

if [ ! -e ${SAMPLE_NAME}.tagcounts.txt ]; then
    
qsub $PROCESS_HOLD -N ".ct${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`

  bash $STAMPIPES/scripts/bwa/tagcounts.bash $SAMPLE_NAME $SAMPLE_NAME.sorted.bam $SAMPLE_NAME.tagcounts.txt
  # upload all data to the LIMS
  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --countsfile ${SAMPLE_NAME}.tagcounts.txt
      
__SCRIPT__

fi

if [ ! -e ${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out -o ! -e ${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt ]; then
    
qsub $PROCESS_HOLD -N ".sp${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`

  make -f $STAMPIPES/makefiles/SPOT/spot-R1.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
  #make -f $STAMPIPES/makefiles/SPOT/spot-Rboth.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
  # upload all data to the LIMS
  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --spotfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out \
      --dupfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt
__SCRIPT__

fi

if [ ! -e ${SAMPLE_NAME}.${GENOME}.bw ]; then

qsub $PROCESS_HOLD -N ".den${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`

  make -f $STAMPIPES/makefiles/densities/density.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
__SCRIPT__

fi
