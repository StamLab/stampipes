# If RETAIN_ORIGINALS is set, then do not delete files we are collating from
# If REDO_COLLATION is set, then do collation again; this will only work if originals still exist

# Ensure that script failures have the script quit before it deletes files
set -e

cd $FASTQ_DIR

FASTQ_NAME=${FLOWCELL}_${SAMPLE_NAME}

echo "Collating $FASTQ_DIR/$FASTQ_NAME"

R1_NUM_FILES=$(find . -maxdepth 1 -name "${SAMPLE_NAME}_R1_???.fastq.gz" | wc -l)

if [[ "$PAIRED" == "True" ]]; then

  R2_NUM_FILES=$(find . -maxdepth 1 -name "${SAMPLE_NAME}_R2_???.fastq.gz" | wc -l)

  if [[ "$R1_NUM_FILES" -ne "$R2_NUM_FILES" ]]; then
    echo "UNEQUAL NUMBER OF FILES FOR $SAMPLE_NAME IN $FASTQ_DIR"
    exit 1
  fi
fi

R1_FILE=${FASTQ_NAME}_R1.fastq.gz
R2_FILE=${FASTQ_NAME}_R2.fastq.gz

function upload {
  if [[ "$SAMPLE_NAME" == LP* ]] ; then
    # Altcode sample, use dedicated script
    python3 "$STAMPIPES/scripts/altcode/upload_fastq.py" --lane "$FLOWCELL_LANE_ID" --r1 "$R1_FILE" --r2 "$R2_FILE"
  else
    # Regular sample, upload old-style
	  UPLOAD_SCRIPT="python3 $STAMPIPES/scripts/lims/upload_data.py --attach_file_contenttype SequencingData.flowcelllane --attach_file_objectid ${FLOWCELL_LANE_ID} --attach_file_type=gzipped-fastq"
	  $UPLOAD_SCRIPT --attach_file_purpose r1-fastq --attach_file "${R1_FILE}"

	  if [ -e "$R2_FILE" ]; then
	    $UPLOAD_SCRIPT --attach_file_purpose r2-fastq --attach_file "${R2_FILE}"
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
  if [ -n "$RETAIN_ORIGINALS" ] ; then
    cmd="cp"
  fi

  "$cmd" "${SAMPLE_NAME}_R1_001.fastq.gz" "$R1_FILE"
  if [ "$PAIRED" == "True" ] ; then
    "$cmd" "${SAMPLE_NAME}_R2_001.fastq.gz" "$R2_FILE"
  fi

  upload
  exit 0
else

  R1_TMP_FILE=$(mktemp)
  R2_TMP_FILE=$(mktemp)

  if [ -e "$R1_FILE" ]; then
    rm "$R1_FILE"
  fi

  if [ -e "$R2_FILE" ]; then
    rm "$R2_FILE"
  fi

  echo "R1: $R1_FILE"
  echo "R2: $R2_FILE"

  for filenum in $(seq -f "%03g" 1 $R1_NUM_FILES)
  do

    echo "Adding ${filenum} to collated files"

    cat ${SAMPLE_NAME}_R1_${filenum}.fastq.gz >> $R1_TMP_FILE

    if [[ "$PAIRED" == "True" ]]; then
      R2_FILE=${FASTQ_NAME}_R2.fastq.gz
      cat ${SAMPLE_NAME}_R2_${filenum}.fastq.gz  >> $R2_TMP_FILE
    fi

  done

  # Ensure we have valid gzipped files
  gzip -t $R1_TMP_FILE

  if [ ! -s $R1_TMP_FILE ]; then
    echo "ERROR: $R1_TMP_FILE is 0 size"
    exit 1
  fi

  if [ "$PAIRED" == "True" ]; then
    gzip -t $R2_TMP_FILE
    if [ ! -s $R2_TMP_FILE ]; then
      echo "ERROR: $R2_TMP_FILE is 0 size"
      exit 1
    fi
  fi

  rsync "$R1_TMP_FILE" "$R1_FILE"
  # Files created in temp directories do not have appropriate
  # permissions; make sure our collated files can be read by
  # anybody
  chmod 644 $R1_FILE
  rm $R1_TMP_FILE
  if [[ "$PAIRED" == "True" ]] ; then
    rsync "$R2_TMP_FILE" "$R2_FILE"
    chmod 644 $R2_FILE
    rm $R2_TMP_FILE
  fi
fi

upload

# Remove existing pre-collation files
if [ ! -n "$RETAIN_ORIGINALS" ]; then
  echo "Removing R1 originals"
  rm ${SAMPLE_NAME}_R1_???.fastq.gz

  if [[ "$PAIRED" == "True" ]]; then
    echo "Removing R2 originals"
    rm ${SAMPLE_NAME}_R2_???.fastq.gz
  fi
fi
