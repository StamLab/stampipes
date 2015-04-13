# Dependencies
source $MODULELOAD
module load bedops/2.4.2
module load bedtools/2.16.2
module load bwa/0.7.0
module load java/jdk1.7.0_05
module load picard/1.118
module load anaconda/2.1.0-2.7
module load samtools/0.1.19
module load gcc/4.7.2
module load R/3.1.0
module load git/2.3.3

FINAL_BAM=${SAMPLE_NAME}.sorted.bam
UNIQUES_BAM=${SAMPLE_NAME}.uniques.sorted.bam

ADAPTER_FILE=${SAMPLE_NAME}.adapters.txt
VERSION_FILE=${SAMPLE_NAME}.versions.txt

bash $STAMPIPES/scripts/versions.bash &> $VERSION_FILE

if [[ ( -n "$ADAPTER_P7" ) && ( -n "ADAPTER_P5" ) ]] ; then
  echo -e "P7\t$ADAPTER_P7\nP5\t$ADAPTER_P5" > $ADAPTER_FILE
fi

if [ ! -e ${SAMPLE_NAME}_R1_fastqc -o ! -e ${SAMPLE_NAME}_R2_fastqc -o -n "$FORCE_FASTQ" ]; then
qsub -N ".fq${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`
  
  make -f $STAMPIPES/makefiles/fastqc.mk
__SCRIPT__
fi

NUMBER_FASTQ_FILES=`find . -maxdepth 1 -name "${SAMPLE_NAME}_R1_???.fastq.gz" | wc -l`
FASTQ_PAIR_HOLDS=""
FASTQ_PAIR_BAMS=""

# Indicate we have started this alignment and upload pertinent information
if [ ! -e ${FINAL_BAM} ]; then

python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
  -t ${LIMS_API_TOKEN} \
  --alignment_id ${ALIGNMENT_ID} \
  --start_alignment_progress \
  --adapter_file $ADAPTER_FILE \
  --version_file $VERSION_FILE

fi

for filenum in $(seq -f "%03g" 1 $NUMBER_FASTQ_FILES)
do
  NAME=".aln${SAMPLE_NAME}_${filenum}_${FLOWCELL}"
  BAMFILE="${SAMPLE_NAME}_${filenum}.sorted.bam"
    
if [ ! -e $BAMFILE -a ! -e ${FINAL_BAM} ]; then

qsub -l h_data=5650M -N ${NAME} -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`
  
  make -f $STAMPIPES/makefiles/bwa/bwa_paired_trimmed.mk \
    FASTQ1_FILE=${SAMPLE_NAME}_R1_${filenum}.fastq.gz \
    FASTQ2_FILE=${SAMPLE_NAME}_R2_${filenum}.fastq.gz \
    OUTBAM=${BAMFILE} \
    TRIMSTATS=${SAMPLE_NAME}_${filenum}.trimstats.txt
__SCRIPT__

  # Only hold on alignments that are being run
  FASTQ_PAIR_HOLDS="$FASTQ_PAIR_HOLDS,$NAME"

fi

# Need to keep track of these even if they have already finshed
# for proper merging
FASTQ_PAIR_BAMS="${BAMFILE} ${FASTQ_PAIR_BAMS}"

done


if [ -n "$FASTQ_PAIR_HOLDS" ]; then
    SPLIT_ALIGN_HOLD="-hold_jid $FASTQ_PAIR_HOLDS"
fi

if [ ! -e ${FINAL_BAM} -o ! -e ${UNIQUES_BAM} ]; then

PROCESS_HOLD="-hold_jid .pb${SAMPLE_NAME}_${FLOWCELL}"

JOBNAME=".pb${SAMPLE_NAME}_${FLOWCELL}"
PROCESSING="$PROCESSING,$JOBNAME"
    
qsub ${SPLIT_ALIGN_HOLD} -N "$JOBNAME" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
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
  make -f $STAMPIPES/makefiles/picard/dups.mk

  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
    -t ${LIMS_API_TOKEN} \
    -f ${FLOWCELL} \
    --alignment_id ${ALIGNMENT_ID} \
    --flowcell_lane_id ${FLOWCELL_LANE_ID} \
    --insertsfile ${SAMPLE_NAME}.CollectInsertSizeMetrics.picard \
    --dupsfile ${SAMPLE_NAME}.MarkDuplicates.picard
  
  if [ "$NUMBER_FASTQ_FILES" -gt "1" ]
  then
    rm $FASTQ_PAIR_BAMS
  fi
__SCRIPT__

fi

if [ ! -e ${SAMPLE_NAME}.tagcounts.txt -o -n "$FORCE_COUNTS" ]; then

JOBNAME=".ct${SAMPLE_NAME}_${FLOWCELL}"
PROCESSING="$PROCESSING,$JOBNAME"
   
qsub $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
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

JOBNAME=".sp${SAMPLE_NAME}_${FLOWCELL}"
PROCESSING="$PROCESSING,$JOBNAME" 
    
qsub $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`

  make -f $STAMPIPES/makefiles/SPOT/spot-R1-paired.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
  # upload all data to the LIMS
  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --spotfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out \
      --spotdupfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt
__SCRIPT__

fi

if [ ! -e ${SAMPLE_NAME}.75_20.${GENOME}.bw ]; then

JOBNAME=".den${SAMPLE_NAME}_${FLOWCELL}"
PROCESSING="$PROCESSING,$JOBNAME" 

qsub $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`

  make -f $STAMPIPES/makefiles/densities/density.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
__SCRIPT__

fi

if [ -e ${PROCESSING} ]; then

qsub -jid_hold ${PROCESSING} -N ".com${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`

  bash $STAMPIPES/scripts/bwa/checkcomplete.bash

  python $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
    -t ${LIMS_API_TOKEN} \
    -f ${FLOWCELL} \
    --alignment_id ${ALIGNMENT_ID} \
    --finish_alignment

__SCRIPT__

fi
