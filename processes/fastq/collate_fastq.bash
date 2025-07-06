# If RETAIN_ORIGINALS is set, then do not delete files we are collating from
# If REDO_COLLATION is set, then do collation again; this will only work if originals still exist

# Ensure that script failures have the script quit before it deletes files
set -e

CLUSTER_NAME=$(scontrol show config | awk '$1 == "ClusterName" {print $3}')
if [[ "$CLUSTER_NAME" == "altius-gene" ]]; then
  module load apptainer/1.3.3
  echo "# Using apptainer"
  export APX="apptainer exec --bind /net/seq/data2/sequencers,/net/seq/data2/flowcells,$STAMPIPES $STAMPIPES/containers/fastq/fastq.sif"
else
  echo "# Not using apptainer"
  export APX=
fi

cd $FASTQ_DIR

INPUT_PREFIX=$SAMPLE_NAME

FASTQ_NAME=${FLOWCELL}_${SAMPLE_NAME}

echo "Collating $FASTQ_DIR/$FASTQ_NAME"

R1_NUM_FILES=$(find . -maxdepth 1 -name "${INPUT_PREFIX}_R1_*.fastq.gz" | wc -l)
R2_NUM_FILES=$(find . -maxdepth 1 -name "${INPUT_PREFIX}_R2_*.fastq.gz" | wc -l)
R3_NUM_FILES=$(find . -maxdepth 1 -name "${INPUT_PREFIX}_R3_*.fastq.gz" | wc -l)
R4_NUM_FILES=$(find . -maxdepth 1 -name "${INPUT_PREFIX}_R4_*.fastq.gz" | wc -l)

R1_FILE=${FASTQ_NAME}_R1.fastq.gz
R2_FILE=${FASTQ_NAME}_R2.fastq.gz
R3_FILE=${FASTQ_NAME}_R3.fastq.gz
R4_FILE=${FASTQ_NAME}_R4.fastq.gz

function upload {
  if [[ "$SAMPLE_NAME" == LP* ]]; then
    # Altcode sample, use dedicated script
    $APX python3 "$STAMPIPES/scripts/altcode/upload_fastq.py" --lane "$FLOWCELL_LANE_ID" --r1 "$R1_FILE" --r2 "$R2_FILE"
  else
    # Regular sample, upload old-style
    UPLOAD_SCRIPT="$APX python3 $STAMPIPES/scripts/lims/upload_data.py --attach_file_contenttype SequencingData.flowcelllane --attach_file_objectid ${FLOWCELL_LANE_ID} --attach_file_type=gzipped-fastq"
    if [ -e "$R1_FILE" ]; then
      $UPLOAD_SCRIPT --attach_file_purpose r1-fastq --attach_file "${R1_FILE}"
    fi

    if [ -e "$R2_FILE" ]; then
      $UPLOAD_SCRIPT --attach_file_purpose r2-fastq --attach_file "${R2_FILE}"
    fi

    if [ -e "$R3_FILE" ]; then
      $UPLOAD_SCRIPT --attach_file_purpose r3-fastq --attach_file "${R3_FILE}"
    fi

    if [ -e "$R4_FILE" ]; then
      $UPLOAD_SCRIPT --attach_file_purpose r4-fastq --attach_file "${R4_FILE}"
    fi
  fi
}

if [ -e "$R1_FILE" -a ! -n "$REDO_COLLATION" ]; then
  echo "Collated files already exist, syncing with LIMS"
  upload
  exit
fi

if [ -n "$REDO_COLLATION" -a "$R1_NUM_FILES" -eq 0 ]; then
  echo "Cannot redo collation, no original files"
  upload
  exit
fi

if [ "$R1_NUM_FILES" -lt "1" ]; then
  echo "No R1 files to work with"
  exit 1
fi

# If only one file, just mv or cp as appropriate
if [ "$R1_NUM_FILES" -eq "1" ]; then
  cmd="mv"
  if [ -n "$RETAIN_ORIGINALS" ]; then
    cmd="cp"
  fi

  if [ "$R1_NUM_FILES" -gt 0 ]; then
    "$cmd" "${INPUT_PREFIX}"_R1*.fastq.gz "$R1_FILE"
  fi
  if [ "$R2_NUM_FILES" -gt 0 ]; then
    "$cmd" "${INPUT_PREFIX}"_R2*.fastq.gz "$R2_FILE"
  fi
  if [ "$R3_NUM_FILES" -gt 0 ]; then
    "$cmd" "${INPUT_PREFIX}"_R3*.fastq.gz "$R3_FILE"
  fi
  if [ "$R4_NUM_FILES" -gt 0 ]; then
    "$cmd" "${INPUT_PREFIX}"_R4*.fastq.gz "$R4_FILE"
  fi

  upload
  exit 0
else

  R1_TMP_FILE=$(mktemp)
  R2_TMP_FILE=$(mktemp)
  R3_TMP_FILE=$(mktemp)
  R4_TMP_FILE=$(mktemp)

  if [ -e "$R1_FILE" ]; then
    rm "$R1_FILE"
  fi

  if [ -e "$R2_FILE" ]; then
    rm "$R2_FILE"
  fi

  if [ -e "$R3_FILE" ]; then
    rm "$R3_FILE"
  fi

  if [ -e "$R4_FILE" ]; then
    rm "$R4_FILE"
  fi

  echo "R1: $R1_FILE"
  echo "R2: $R2_FILE"
  echo "R3: $R3_FILE"
  echo "R4: $R4_FILE"

  for filenum in $(seq -f "%03g" 1 $R1_NUM_FILES); do

    echo "Adding ${filenum} to collated files"

    if [ -e "${SAMPLE_NAME}_R1_${filenum}.fastq.gz" ]; then
      cat ${SAMPLE_NAME}_R1_${filenum}.fastq.gz >>$R1_TMP_FILE
    fi
    if [ -e "${SAMPLE_NAME}_R2_${filenum}.fastq.gz" ]; then
      cat ${SAMPLE_NAME}_R2_${filenum}.fastq.gz >>$R2_TMP_FILE
    fi
    if [ -e "${SAMPLE_NAME}_R3_${filenum}.fastq.gz" ]; then
      cat ${SAMPLE_NAME}_R3_${filenum}.fastq.gz >>$R3_TMP_FILE
    fi
    if [ -e "${SAMPLE_NAME}_R4_${filenum}.fastq.gz" ]; then
      cat ${SAMPLE_NAME}_R4_${filenum}.fastq.gz >>$R4_TMP_FILE
    fi

  done

  # Ensure we have valid gzipped files
  if [ -s $R1_TMP_FILE ]; then
    gzip -t $R1_TMP_FILE
  else
    echo "ERROR: $R1_TMP_FILE is 0 size"
    exit 1
  fi

  if [ -s $R2_TMP_FILE ]; then
    gzip -t $R2_TMP_FILE
  elif [ "$R2_NUM_FILES" -gt 0 ]; then
    echo "ERROR: $R2_TMP_FILE is 0 size"
    exit 1
  fi

  if [ -s $R3_TMP_FILE ]; then
    gzip -t $R3_TMP_FILE
  elif [ "$R3_NUM_FILES" -gt 0 ]; then
    echo "ERROR: $R3_TMP_FILE is 0 size"
    exit 1
  fi

  if [ -s $R4_TMP_FILE ]; then
    gzip -t $R4_TMP_FILE
  elif [ "$R4_NUM_FILES" -gt 0 ]; then
    echo "ERROR: $R4_TMP_FILE is 0 size"
    exit 1
  fi

  # Move temp files to final locations with proper permissions
  if [ -s $R1_TMP_FILE ]; then
    rsync "$R1_TMP_FILE" "$R1_FILE"
    chmod 644 $R1_FILE
  fi
  rm $R1_TMP_FILE

  if [ -s $R2_TMP_FILE ]; then
    rsync "$R2_TMP_FILE" "$R2_FILE"
    chmod 644 $R2_FILE
  fi
  rm $R2_TMP_FILE

  if [ -s $R3_TMP_FILE ]; then
    rsync "$R3_TMP_FILE" "$R3_FILE"
    chmod 644 $R3_FILE
  fi
  rm $R3_TMP_FILE

  if [ -s $R4_TMP_FILE ]; then
    rsync "$R4_TMP_FILE" "$R4_FILE"
    chmod 644 $R4_FILE
  fi
  rm $R4_TMP_FILE
fi

upload

# Remove existing pre-collation files
if [ ! -n "$RETAIN_ORIGINALS" ]; then
  if [ "$R1_NUM_FILES" -gt 0 ]; then
    echo "Removing R1 originals"
    rm ${SAMPLE_NAME}_R1_???.fastq.gz
  fi

  if [ "$R2_NUM_FILES" -gt 0 ]; then
    echo "Removing R2 originals"
    rm ${SAMPLE_NAME}_R2_???.fastq.gz
  fi

  if [ "$R3_NUM_FILES" -gt 0 ]; then
    echo "Removing R3 originals"
    rm ${SAMPLE_NAME}_R3_???.fastq.gz
  fi

  if [ "$R4_NUM_FILES" -gt 0 ]; then
    echo "Removing R4 originals"
    rm ${SAMPLE_NAME}_R4_???.fastq.gz
  fi
fi
