# Unpaired bwa alignment, untrimmed

# Dependencies
source $MODULELOAD
module load bedops/2.4.19
module load bedtools/2.25.0
module load bwa/0.7.12
module load jdk/1.8.0_92
module load picard/1.120
module load samtools/1.3
module load gcc/4.7.2
module load R/3.2.5
module load git/2.3.3
module load coreutils/8.25
module load pigz/2.3.3

# Load in this order specifically, currently the python3 activation
# overwrites the default "python" call, against advice
module load python/3.5.1
module load pysam/0.9.0
module load python/2.7.11

MAX_MISMATCHES=2
MIN_MAPPING_QUALITY=10
MAX_INSERT_SIZE=750

JOB_BASENAME=${SAMPLE_NAME}_${FLOWCELL}_ALIGN#${ALIGNMENT_ID}

export FINAL_BAM=${SAMPLE_NAME}.sorted.bam
export UNIQUES_BAM=${SAMPLE_NAME}.uniques.sorted.bam

export VERSION_FILE=${SAMPLE_NAME}.versions.txt

export FASTQ_TMP=$ALIGN_DIR/fastq

cd $ALIGN_DIR

if [ -n "$REDO_ALIGNMENT" ]; then
    bash $STAMPIPES/scripts/bwa/reset_alignment.bash
fi

bash $STAMPIPES/scripts/versions.bash &> $VERSION_FILE

# Indicate we have started this alignment and upload pertinent information
if [ ! -e ${FINAL_BAM} ]; then

if [ -n ${ALIGNMENT_ID} ]; then
python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
  -t ${LIMS_API_TOKEN} \
  --alignment_id ${ALIGNMENT_ID} \
  --start_alignment_progress \
  --version_file $VERSION_FILE
fi

if [ ! -e "$FASTQ_TMP/${SAMPLE_NAME}_R1_000.fastq.gz" ]; then
  bash $STAMPIPES/scripts/fastq/splitfastq.bash $FASTQ_TMP $R1_FASTQ $R2_FASTQ
fi

fi

NUMBER_FASTQ_FILES=`find $FASTQ_TMP -maxdepth 1 -name "${SAMPLE_NAME}_R1_???.fastq.gz" | wc -l`
FASTQ_PAIR_HOLDS=""
FASTQ_PAIR_BAMS=""

for filenum in $(seq -f "%03g" 0 `expr $NUMBER_FASTQ_FILES - 1`)
do
  NAME=".aln${JOB_BASENAME}_${filenum}"
  BAMFILE="${SAMPLE_NAME}_${filenum}.sorted.bam"
    
if [ ! -e $BAMFILE -a ! -e ${FINAL_BAM} ]; then

qsub -N ${NAME} -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  make -f $STAMPIPES/makefiles/bwa/bwa_single.mk \
    FASTQ_FILE=${FASTQ_TMP}/${SAMPLE_NAME}_R1_${filenum}.fastq.gz \
    OUTBAM=${BAMFILE}

  echo "FINISH: " `date`

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

if [ ! -e ${FINAL_BAM}.bai -o ! -e ${UNIQUES_BAM}.bai ]; then

# If we are redoing this part, then we should make sure
# to redo all the other results as well
export FORCE_COUNTS=1
export REPROCESS=1

PROCESS_HOLD="-hold_jid .pb${JOB_BASENAME}"

JOBNAME=".pb${JOB_BASENAME}"
PROCESSING="$PROCESSING,$JOBNAME"

qsub ${SPLIT_ALIGN_HOLD} ${SUBMIT_SLOTS} -N "$JOBNAME" -V -cwd -S /bin/bash << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  if [ ! -e ${FINAL_BAM} ]; then
    echo "START MERGE: "
    date
  
    if [ "$NUMBER_FASTQ_FILES" -eq "1" ]; then
      mv ${SAMPLE_NAME}_000.sorted.bam ${FINAL_BAM}
    else
      samtools merge ${FINAL_BAM} ${FASTQ_PAIR_BAMS}
    fi
  
    echo "FINISH MERGE: "
    date
  else
    echo "$FINAL_BAM exists already"
  fi
  
  echo "START PROCESS BAM: "
  date

  make -f $STAMPIPES/makefiles/bwa/process_unpaired_bam.mk

  echo "FINISH PROCESS BAM: "
  date

  echo "START PICARD DUP: "
  date

  make -f $STAMPIPES/makefiles/picard/dups.mk

  echo "FINISH PICARD DUP: " date

  python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
    -t ${LIMS_API_TOKEN} \
    -f ${FLOWCELL} \
    --alignment_id ${ALIGNMENT_ID} \
    --flowcell_lane_id ${FLOWCELL_LANE_ID} \
    --dupsfile ${SAMPLE_NAME}.MarkDuplicates.picard
  
  if [ "$NUMBER_FASTQ_FILES" -gt "1" ]
  then
    rm -f $FASTQ_PAIR_BAMS
  fi

  echo "FINISH: "
  date

__SCRIPT__

fi

if [ ! -e ${SAMPLE_NAME}.tagcounts.txt -o -n "$FORCE_COUNTS" ]; then

  JOBNAME=".ct${JOB_BASENAME}"
  PROCESSING="$PROCESSING,$JOBNAME"

  qsub $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: "
  hostname

  echo "START: " 
  date

  time bash $STAMPIPES/scripts/bwa/tagcounts.bash $SAMPLE_NAME $SAMPLE_NAME.sorted.bam $SAMPLE_NAME.tagcounts.txt $FASTQ_DIR
  # upload all data to the LIMS
  python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --countsfile ${SAMPLE_NAME}.tagcounts.txt

  echo "FINISH: "
  date
      
__SCRIPT__

fi

if [ ! -e ${SAMPLE_NAME}.rand.uniques.sorted.spot.out -o ! -e ${SAMPLE_NAME}.rand.uniques.sorted.spotdups.txt ]; then

JOBNAME=".sp${JOB_BASENAME}"
PROCESSING="$PROCESSING,$JOBNAME" 
    
qsub $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  make -f $STAMPIPES/makefiles/SPOT/spot-single.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
  # upload all data to the LIMS
  python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
      -t ${LIMS_API_TOKEN} \
      -f ${FLOWCELL} \
      --alignment_id ${ALIGNMENT_ID} \
      --flowcell_lane_id ${FLOWCELL_LANE_ID} \
      --spotfile ${SAMPLE_NAME}.rand.uniques.sorted.spot.out \
      --spotdupfile ${SAMPLE_NAME}.rand.uniques.sorted.spotdups.txt

  echo "FINISH: "
  date

__SCRIPT__

fi

if [ ! -e ${SAMPLE_NAME}.75_20.${GENOME}.bw ]; then

  JOBNAME=".den${JOB_BASENAME}"
  PROCESSING="$PROCESSING,$JOBNAME"

  qsub $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  make -f $STAMPIPES/makefiles/densities/density.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
    READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME

  echo "FINISH: "
  date

__SCRIPT__

fi

if [ -n "$PROCESSING" ]; then

qsub -hold_jid ${PROCESSING} -N ".com${JOB_BASENAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: "
  hostname
  
  echo "START: "
  date

  bash $STAMPIPES/scripts/bwa/unpairedcheckcomplete.bash

  python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
    -t ${LIMS_API_TOKEN} \
    -f ${FLOWCELL} \
    --alignment_id ${ALIGNMENT_ID} \
    --finish_alignment

  bash $STAMPIPES/scripts/bwa/attachfiles.bash

  rm -rf $ALIGN_DIR/fastq

  echo "FINISH: "
  date
__SCRIPT__

fi
