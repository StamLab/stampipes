# Dependencies
source "$MODULELOAD"
module load bedops/2.4.19
module load bedtools/2.25.0
module load bwa/0.7.12
module load jdk/1.8.0_92
module load picard/1.120
module load samtools/1.3
module load gcc/4.7.2
module load R/3.2.5

module load python/3.5.1
module load pysam/0.9.0

FINAL_BAM=${SAMPLE_NAME}.sorted.bam
UNIQUES_BAM=${SAMPLE_NAME}.uniques.sorted.bam

export FASTQ_TMP=$ALIGN_DIR/fastq

set -x

bash "$STAMPIPES/scripts/versions.bash" &>"$SAMPLE_NAME.versions.txt"

if [[ ! -e "$FINAL_BAM" ]]; then
  if [ ! -e "$FASTQ_TMP/${SAMPLE_NAME}_R1_000.fastq.gz" ]; then
    bash "$STAMPIPES/scripts/fastq/splitfastq.bash" "$FASTQ_TMP" "$R1_FASTQ" "$R2_FASTQ"
  fi

  FASTQ_PAIR_HOLDS=""
  FASTQ_PAIR_BAMS=""
  NUMBER_FASTQ_FILES=$(find "$FASTQ_TMP" -maxdepth 1 -name "${SAMPLE_NAME}_R1_???.fastq.gz" | wc -l)
  for filenum in $(seq -f "%03g" 0 $((NUMBER_FASTQ_FILES - 1))); do
    NAME=".aln${SAMPLE_NAME}_${filenum}_${FLOWCELL}"
    BAMFILE="${SAMPLE_NAME}_${filenum}.sorted.bam"

    if [[ ! -e "$BAMFILE" ]]; then
      qsub -N "$NAME" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
      set -x -e -o pipefail
      echo "Hostname: " \$(hostname)

      make -f "$STAMPIPES/makefiles/bwa/bwa_paired.mk" \
        FASTQ1_FILE="$FASTQ_TMP/${SAMPLE_NAME}_R1_${filenum}.fastq.gz" \
        FASTQ2_FILE="$FASTQ_TMP/${SAMPLE_NAME}_R2_${filenum}.fastq.gz" \
        OUTBAM="$BAMFILE" \
        ADAPTERFILE="$STAMPIPES_DATA/adapters/default.adapters"
__SCRIPT__

      # Only hold on alignments that are being run
      FASTQ_PAIR_HOLDS="$FASTQ_PAIR_HOLDS,$NAME"

      # Need to keep track of these even if they have already finshed
      # for proper merging
      FASTQ_PAIR_BAMS="${BAMFILE} ${FASTQ_PAIR_BAMS}"

    fi

  done
fi

if [ -n "$FASTQ_PAIR_HOLDS" ]; then
  HOLD="-hold_jid $FASTQ_PAIR_HOLDS"
fi

if [[ ! -e "$FINAL_BAM.bai" || ! -e "$UNIQUES_BAM.bai" ]]; then

  PROCESS_HOLD="-hold_jid .pb${SAMPLE_NAME}_${FLOWCELL}"

  qsub $HOLD -N ".pb${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " \$(hostname)

  if [ ! -e "$FINAL_BAM" ]; then
    if [ "$NUMBER_FASTQ_FILES" -eq "1" ]
    then
      mv "${SAMPLE_NAME}_000.sorted.bam" "${SAMPLE_NAME}.sorted.bam"
    else
      samtools merge "${FINAL_BAM}" ${FASTQ_PAIR_BAMS}
    fi
  fi

  make -f $STAMPIPES/makefiles/bwa/process_paired_bam.mk
  make -f $STAMPIPES/makefiles/picard/dups.mk

  python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
    -t ${LIMS_API_TOKEN} \
    -f ${FLOWCELL} \
    --alignment_id ${ALIGNMENT_ID} \
    --flowcell_lane_id ${FLOWCELL_LANE_ID} \
    --insertsfile ${SAMPLE_NAME}.CollectInsertSizeMetrics.picard \
    --dupsfile ${SAMPLE_NAME}.MarkDuplicates.picard

  if [ "$NUMBER_FASTQ_FILES" -gt "1" ]
  then
    rm -f $FASTQ_PAIR_BAMS
  fi
__SCRIPT__

fi

if [[ ! -e "${SAMPLE_NAME}.tagcounts.txt" || -n "$FORCE_COUNTS" ]]; then

  qsub $PROCESS_HOLD -N ".ct${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " \$(hostname)

  bash "$STAMPIPES/scripts/bwa/tagcounts.bash" "$SAMPLE_NAME" "${SAMPLE_NAME}.sorted.bam" "${SAMPLE_NAME}.tagcounts.txt"
  # upload all data to the LIMS
  python3 $STAMPIPES/scripts/lims/upload_data.py -a "$LIMS_API_URL" \
    -t "${LIMS_API_TOKEN}" \
    -f "${FLOWCELL}" \
    --alignment_id "${ALIGNMENT_ID}" \
    --flowcell_lane_id "${FLOWCELL_LANE_ID}" \
    --countsfile "${SAMPLE_NAME}.tagcounts.txt"

__SCRIPT__

fi

if [[ ! -e "${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out" || ! -e "${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt" ]]; then

  qsub $PROCESS_HOLD -N ".sp${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " \$(hostname)

  # SPOT process requires python 2
  source "$MODULELOAD"
  module load python/2.7.11

  make -f"" $STAMPIPES/makefiles/SPOT/spot-R1-paired.mk "BWAINDEX=$BWAINDEX" "ASSAY=$ASSAY" "GENOME=$GENOME" \
    "READLENGTH=$READLENGTH" "SAMPLE_NAME=$SAMPLE_NAME"
  # upload all data to the LIMS
  python3 "$STAMPIPES/scripts/lims/upload_data.py" -a "$LIMS_API_URL" \
    -t "$LIMS_API_TOKEN" \
    -f "$FLOWCELL" \
    --alignment_id "$ALIGNMENT_ID" \
    --flowcell_lane_id "$FLOWCELL_LANE_ID" \
    --spotfile "${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out" \
    --spotdupfile "${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt"
__SCRIPT__

fi

if [ ! -e "${SAMPLE_NAME}.75_20.${GENOME}.bw" ]; then

  qsub $PROCESS_HOLD -N ".den${SAMPLE_NAME}_${FLOWCELL}" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " \$(hostname)

  make -f "$STAMPIPES/makefiles/densities/density.mk" "BWAINDEX=$BWAINDEX" "ASSAY=$ASSAY" "GENOME=$GENOME" \
    "READLENGTH=$READLENGTH" "SAMPLE_NAME=$SAMPLE_NAME"
__SCRIPT__

fi
