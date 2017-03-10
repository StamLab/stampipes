#!/bin/bash

# Dependencies
source "$MODULELOAD"
module load zlib/1.2.8
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

module load modwt/1.0
module load hotspot2/20161006

# Load in this order specifically, currently the python3 activation
# overwrites the default "python" call, against advice
module load python/3.5.1
module load pysam/0.9.0
module load python/2.7.11

export MAX_MISMATCHES=2
export MIN_MAPPING_QUALITY=10
export MAX_INSERT_SIZE=750

export BASE_PRIORITY=${BASE_PRIORITY:-"-10"}

if [[ -z "$ALN_PRIORITY" ]]; then
  export ALN_PRIORITY=$((BASE_PRIORITY - 10))
fi

JOB_BASENAME=${SAMPLE_NAME}_${FLOWCELL}_ALIGN#${ALIGNMENT_ID}

export FINAL_BAM=${SAMPLE_NAME}.sorted.bam
export UNIQUES_BAM=${SAMPLE_NAME}.uniques.sorted.bam

export ADAPTER_FILE=${SAMPLE_NAME}.adapters.txt
export VERSION_FILE=${SAMPLE_NAME}.versions.txt

export FASTQ_TMP=$ALIGN_DIR/fastq


cd "$ALIGN_DIR"

if [[ -n "$REDO_ALIGNMENT" ]]; then
  bash "$STAMPIPES/scripts/bwa/reset_alignment.bash"
fi

bash "$STAMPIPES/scripts/versions.bash" &>"$VERSION_FILE"

if [[ (-n "$ADAPTER_P7") && (-n "ADAPTER_P5") ]]; then
  echo -e "P7\t$ADAPTER_P7\nP5\t$ADAPTER_P5" >"$ADAPTER_FILE"
fi

if [[ "$TRIM_READS_TO" -gt 0 && "$READLENGTH" -gt "$TRIM_READS_TO" ]]; then
  NEED_TRIMMING=TRUE
  READLENGTH=$TRIM_READS_TO
fi

# Indicate we have started this alignment and upload pertinent information
if [[ ! -e "$FINAL_BAM" ]]; then

  python3 "$STAMPIPES/scripts/lims/upload_data.py" -a "$LIMS_API_URL" \
    -t "$LIMS_API_TOKEN" \
    --alignment_id "$ALIGNMENT_ID" \
    --start_alignment_progress \
    --adapter_file "$ADAPTER_FILE" \
    --version_file "$VERSION_FILE"

  if [[ ! -e "$FASTQ_TMP/${SAMPLE_NAME}_R1_000.fastq.gz" ]]; then
    bash "$STAMPIPES/scripts/fastq/splitfastq.bash" "$FASTQ_TMP" "$R1_FASTQ" "$R2_FASTQ"
  fi

fi

mkdir -p "$FASTQ_TMP"
NUMBER_FASTQ_FILES=$(find "$FASTQ_TMP" -maxdepth 1 -name "${SAMPLE_NAME}_R1_???.fastq.gz" | wc -l)
FASTQ_PAIR_HOLDS=""
FASTQ_PAIR_BAMS=""

for filenum in $(seq -f "%03g" 0 $((NUMBER_FASTQ_FILES - 1))); do
  NAME=".aln${JOB_BASENAME}_${filenum}"
  BAMFILE="${SAMPLE_NAME}_${filenum}.sorted.bam"

  if [[ ! -e "$BAMFILE" && ! -e "$FINAL_BAM" ]]; then

    qsub -p "$ALN_PRIORITY" -N "$NAME" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
      set -x -e -o pipefail

      echo "Hostname: "
      hostname

      echo "START: "
      date

      fastq1=${FASTQ_TMP}/${SAMPLE_NAME}_R1_${filenum}.fastq.gz
      fastq2=${FASTQ_TMP}/${SAMPLE_NAME}_R2_${filenum}.fastq.gz

      if [[ -n "$NEED_TRIMMING" ]] ; then

        # Define trimming function
        trim_to_length(){
          echo "Trimming \$1 to \$2 bp..."
          mv "\$1" "\$1.tmp"
          zcat \$1.tmp \
            | awk 'NR%2==0 {print substr(\$0, 1, $TRIM_TO)} NR%2!=0' \
            | gzip -c \
            > \$1
          rm \$1.tmp
        }

        trim_to_length \$fastq1 "$TRIM_READS_TO"
        trim_to_length \$fastq2 "$TRIM_READS_TO"
      fi

      # If we have UMI info, create new FASTQ files with the UMI info
      # and use those as input to the makefile
      if [[ -n "$UMI" ]]; then
        new1=\$(mktemp)
        new2=\$(mktemp)

        case "$UMI_METHOD" in
        thruplex)
          time python3 $STAMPIPES/scripts/umi/extract_umt.py \
            <(zcat \$fastq1) \
            <(zcat \$fastq2) \
            >(gzip -c -1 > \$new1) \
            >(gzip -c -1 > \$new2)
          ;;
        *)
          time python3 $STAMPIPES/scripts/umi/fastq_umi_add.py \$fastq1 \$new1
          time python3 $STAMPIPES/scripts/umi/fastq_umi_add.py \$fastq2 \$new2
          ;;
        esac

        fastq1=\$new1
        fastq2=\$new2
      fi

      make -f "$STAMPIPES/makefiles/bwa/bwa_paired_trimmed.mk" \
        "FASTQ1_FILE=\$fastq1" \
        "FASTQ2_FILE=\$fastq2" \
        "OUTBAM=$BAMFILE" \
        "TRIMSTATS=${SAMPLE_NAME}_${filenum}.trimstats.txt"

      echo "FINISH: "
      date

__SCRIPT__

    # Only hold on alignments that are being run
    FASTQ_PAIR_HOLDS="$FASTQ_PAIR_HOLDS,$NAME"

  fi

  # Need to keep track of these even if they have already finshed
  # for proper merging
  FASTQ_PAIR_BAMS="${BAMFILE} ${FASTQ_PAIR_BAMS}"

done

if [[ -n "$FASTQ_PAIR_HOLDS" ]]; then
  SPLIT_ALIGN_HOLD="-hold_jid $FASTQ_PAIR_HOLDS"
fi

if [[ ! -e "$FINAL_BAM.bai" || ! -e "$UNIQUES_BAM.bai" ]]; then

  # If we are redoing this part, then we should make sure
  # to redo all the other results as well
  export FORCE_COUNTS=1
  export REPROCESS=1

  PROCESS_HOLD="-hold_jid .pb${JOB_BASENAME}"

  JOBNAME=".pb${JOB_BASENAME}"
  PROCESSING="$PROCESSING,$JOBNAME"

  # If we are processing UMI, we will need a lot of power for sorting!
  if [[ "$UMI_METHOD" = "single-strand" ]]; then
    export SUBMIT_SLOTS="-pe threads 1"
  fi

  export PRE_DUP_BAM=${SAMPLE_NAME}.predup.sorted.bam
  export PRE_DUP_BAM_PREFIX=${PRE_DUP_BAM%.*}
  export FINAL_BAM_PREFIX=${FINAL_BAM%.*}

  qsub -p $BASE_PRIORITY ${SPLIT_ALIGN_HOLD} ${SUBMIT_SLOTS} -N "$JOBNAME" -V -cwd -S /bin/bash <<__SCRIPT__
    set -x -e -o pipefail

    echo "Hostname: "
    hostname

    echo "START: "
    date

    if [[ ! -e ${FINAL_BAM} ]]; then
      echo "START MERGE: "
      date

      if [[ "$NUMBER_FASTQ_FILES" -eq "1" ]]; then
        mv ${SAMPLE_NAME}_000.sorted.bam ${FINAL_BAM}
      else
        samtools merge ${FINAL_BAM} ${FASTQ_PAIR_BAMS}
      fi

  # Skip Dup marking

      # If we are working with a UMI, we can create the duplicate score and replace
      # the final BAM with one marked with duplicates
      if [[ "$UMI" = "True" ]]; then
        echo "START UMI DUP: "
        date

        rsync $FINAL_BAM \$TMPDIR/${PRE_DUP_BAM}
        rm -f $FINAL_BAM

        make -f $STAMPIPES/makefiles/umi/mark_duplicates.mk INPUT_BAM_FILE=\$TMPDIR/${PRE_DUP_BAM} \
          OUTPUT_BAM_FILE=\$TMPDIR/${FINAL_BAM}.presort
        samtools sort \$TMPDIR/${FINAL_BAM}.presort > $FINAL_BAM_PREFIX.bam

        echo "FINISH UMI DUP: "
        date
      fi

      echo "FINISH MERGE: "
      date
    else
      echo "$FINAL_BAM exists already"
    fi

  # Skip dup marking and inserts and uniquely-mapping steps.
    echo "START PROCESS BAM: "
    date

    make -f $STAMPIPES/makefiles/bwa/process_paired_bam.mk

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
      --insertsfile ${SAMPLE_NAME}.CollectInsertSizeMetrics.picard \
      --dupsfile ${SAMPLE_NAME}.MarkDuplicates.picard

    if [[ "$NUMBER_FASTQ_FILES" -gt "1" ]]
    then
      rm -f $FASTQ_PAIR_BAMS
    fi

    echo "FINISH: "
    date

__SCRIPT__

fi

# Preseq duplicate estimation
if [[ ! -s "$SAMPLE_NAME.uniques.preseq.targets.txt" || ! -s "$SAMPLE_NAME.uniques.dup.lorentz.txt" ]]; then
  JOBNAME=".ps${JOB_BASENAME}"
  PROCESSING="$PROCESSING,$JOBNAME"
  qsub -p "$BASE_PRIORITY" $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash >/dev/stderr <<'__SCRIPT__'
    set -x
    hist=$SAMPLE_NAME.uniques.duphist.txt
    preseq=$SAMPLE_NAME.uniques.preseq.txt
    targets=$SAMPLE_NAME.uniques.preseq.targets.txt
    lorentz=$SAMPLE_NAME.uniques.dup.lorentz.txt

    FRAGMENT_TARGETS=(
      20000000
      30000000
      40000000
     100000000
     150000000
     250000000
  )

    python3 "$STAMPIPES/scripts/bam/mark_dups.py" -i "$UNIQUES_BAM" -o /dev/null --hist "$hist"

    # Get gini / robinhood metrics
    python3 "$STAMPIPES/scripts/utility/lorentz.py" < "$hist" > "$lorentz"
    python3 "$STAMPIPES/scripts/lims/upload_data.py" -a "$LIMS_API_URL" -t "$LIMS_API_TOKEN" --alignment_id "$ALIGNMENT_ID" --countsfile "$lorentz"

    # Get preseq metrics
    preseq lc_extrap -hist "$hist" -extrap 1.001e9 -s 1e6 -v > "$preseq"

    python3 "$STAMPIPES/scripts/lims/upload_data.py" -a "$LIMS_API_URL" -t "$LIMS_API_TOKEN" \
      --attach_file "$preseq" \
      --attach_file_purpose preseq-estimation \
      --attach_file_type plaintext \
      --attach_file_objectid "$ALIGNMENT_ID" \
      --attach_file_contenttype SequencingData.flowcelllanealignment


    # Write preseq targets out
    rm -f "$targets"
    for target in "${FRAGMENT_TARGETS[@]}" ; do
      fragments_needed=$(awk -v "target=$target" 'NR>1 && $2 > target {printf "%d", $1; exit}' "$preseq")
      if [ -n "$fragments_needed" ]; then
        echo -e "preseq-est-for-$target\t$fragments_needed" >> "$targets"
      fi
    done
    maximum=$(awk 'END{print int($2)}' "$preseq")
    echo -e "preseq-est-max\t$maximum" >> "$targets"
    python3 "$STAMPIPES/scripts/lims/upload_data.py" -a "$LIMS_API_URL" -t "$LIMS_API_TOKEN" --alignment_id "$ALIGNMENT_ID" --countsfile "$targets"
__SCRIPT__
fi

if [[ ! -e "$SAMPLE_NAME.uniques.spot2.out" ]]; then
  JOBNAME=".sp2$JOB_BASENAME"
  PROCESSING="$PROCESSING,$JOBNAME"

  HOTSPOT_PREFIX=$(basename "$UNIQUES_BAM" .bam)
  qsub -p "$BASE_PRIORITY" $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
    set -e -u -o errexit
    num_fragments=\$(samtools idxstats "$UNIQUES_BAM" | awk '{x+=\$3}END{print x/2}')
    if [[ \$num_fragments -lt 10000000 ]]; then
      echo "Not enough fragments to get a reliable SPOT2 score, skipping"
      exit
    fi

    hotspot2.sh \
      -c "$BWAINDEX.chrom_sizes.bed" \
      -C "$BWAINDEX.K$READLENGTH.center_sites.n100.starch" \
      -M "$BWAINDEX.K$READLENGTH.mappable_only.bed" \
      -s 12345 \
      -F 0.1 \
      -f 0.05 \
      "$UNIQUES_BAM" \
      "$ALIGN_DIR/hotspot2"

      cd "$ALIGN_DIR"
      cp "hotspot2/$HOTSPOT_PREFIX.SPOT.txt" "$SAMPLE_NAME.uniques.spot2.out"
      "$STAMPIPES/scripts/SPOT/info.sh" \
        "hotspot2/$HOTSPOT_PREFIX.hotspots.fdr0.05.starch" \
        hotspot2 \
        hotspot2/$HOTSPOT_PREFIX.SPOT.txt > "$HOTSPOT_PREFIX.hotspot2.info"

__SCRIPT__

fi

# Skip a bunch of stuff
if [[ ! -e "$SAMPLE_NAME.tagcounts.txt" || -n "$FORCE_COUNTS" ]]; then

  JOBNAME=".ct${JOB_BASENAME}"
  PROCESSING="$PROCESSING,$JOBNAME"

  qsub -p "$BASE_PRIORITY" $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
    set -x -e -o pipefail
    echo "Hostname: "
    hostname

    echo "START: "
    date

    time bash $STAMPIPES/scripts/bwa/tagcounts.bash $SAMPLE_NAME $SAMPLE_NAME.sorted.bam $SAMPLE_NAME.tagcounts.txt $R1_FASTQ $R2_FASTQ
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

if [[ ! -e "$SAMPLE_NAME.R1.rand.uniques.sorted.spot.out" || ! -e "$SAMPLE_NAME.R1.rand.uniques.sorted.spotdups.txt" ]]; then

  JOBNAME=".sp${JOB_BASENAME}"
  PROCESSING="$PROCESSING,$JOBNAME"

  qsub -p $BASE_PRIORITY $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
    set -x -e -o pipefail

    echo "Hostname: "
    hostname

    echo "START: "
    date

    make -f $STAMPIPES/makefiles/SPOT/spot-R1-paired.mk BWAINDEX=$BWAINDEX ASSAY=$ASSAY GENOME=$GENOME \
      READLENGTH=$READLENGTH SAMPLE_NAME=$SAMPLE_NAME
    # upload all data to the LIMS
    python3 $STAMPIPES/scripts/lims/upload_data.py -a ${LIMS_API_URL} \
        -t ${LIMS_API_TOKEN} \
        -f ${FLOWCELL} \
        --alignment_id ${ALIGNMENT_ID} \
        --flowcell_lane_id ${FLOWCELL_LANE_ID} \
        --spotfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spot.out \
        --spotdupfile ${SAMPLE_NAME}.R1.rand.uniques.sorted.spotdups.txt

    echo "FINISH: "
    date

__SCRIPT__

fi

if [[ ! -e "$SAMPLE_NAME.75_20.$GENOME.bw" ]]; then

  JOBNAME=".den${JOB_BASENAME}"
  PROCESSING="$PROCESSING,$JOBNAME"

  qsub -p "$BASE_PRIORITY" $PROCESS_HOLD -N "$JOBNAME" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
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

if [[ -n "$PROCESSING" ]]; then

  qsub -p "$BASE_PRIORITY" -hold_jid "$PROCESSING" -N ".com$JOB_BASENAME" -V -cwd -S /bin/bash >/dev/stderr <<__SCRIPT__
    set -x -e -o pipefail
    echo "Hostname: "
    hostname

    echo "START: "
    date

    bash $STAMPIPES/scripts/bwa/checkcomplete.bash

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
