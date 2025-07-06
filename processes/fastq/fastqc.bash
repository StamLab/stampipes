# Dependencies
[[ -s "$MODULELOAD" ]] && source "$MODULELOAD"
{
  module load jdk/1.8.0_92
  module load picard/2.8.1
  module load fastqc/0.11.5
} || true # ignore module load failures

export FASTQ_NAME=${FLOWCELL}_${SAMPLE_NAME}

export R1_FASTQ=${FASTQ_NAME}_R1.fastq.gz
export R2_FASTQ=${FASTQ_NAME}_R2.fastq.gz
export R3_FASTQ=${FASTQ_NAME}_R3.fastq.gz
export R4_FASTQ=${FASTQ_NAME}_R4.fastq.gz

export R1_FASTQC=${FASTQ_NAME}_R1_fastqc.zip
export R2_FASTQC=${FASTQ_NAME}_R2_fastqc.zip
export R3_FASTQC=${FASTQ_NAME}_R3_fastqc.zip
export R4_FASTQC=${FASTQ_NAME}_R4_fastqc.zip

export TOP_UMIS=${SAMPLE_NAME}.topumis.txt.gz

cd $FASTQ_DIR

CLUSTER_NAME=$(scontrol show config | awk '$1 == "ClusterName" {print $3}')
if [[ "$CLUSTER_NAME" == "altius-gene" ]]; then
  module load apptainer/1.3.3
  echo "# Using apptainer"
  export APX="apptainer exec --bind /net/seq/data2/sequencers,/net/seq/data2/flowcells,$STAMPIPES $STAMPIPES/containers/fastq/fastq.sif"
else
  echo "# Not using apptainer"
  export APX=
fi

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: "
date

cd $FASTQ_DIR
if [ -e "$R1_FASTQ" ]; then
  $APX make -f $STAMPIPES/makefiles/fastqc.mk FASTQ_FILE=$R1_FASTQ FASTQC_FILE=$R1_FASTQC
fi
if [ -e "$R2_FASTQ" ]; then
  $APX make -f $STAMPIPES/makefiles/fastqc.mk FASTQ_FILE=$R2_FASTQ FASTQC_FILE=$R2_FASTQC
fi
if [ -e "$R3_FASTQ" ]; then
  $APX make -f $STAMPIPES/makefiles/fastqc.mk FASTQ_FILE=$R3_FASTQ FASTQC_FILE=$R3_FASTQC
fi
if [ -e "$R4_FASTQ" ]; then
  $APX make -f $STAMPIPES/makefiles/fastqc.mk FASTQ_FILE=$R4_FASTQ FASTQC_FILE=$R4_FASTQC
fi

if [ "$UMI" = "True" ]; then
  echo "Tallying up top UMI tags seen in R1"
  zcat ${R1_FASTQ} | grep "^@" | cut -f 2 -d "+" | sort | uniq -c | sort -n -r | gzip -c >${TOP_UMIS}
fi

# Build upload command dynamically based on which FastQC files exist
UPLOAD_CMD="$APX python3 ${STAMPIPES}/scripts/lims/upload_data_withr3.py -f ${FLOWCELL} --flowcell_lane_id=${FLOWCELL_LANE_ID}"
if [ -e "$R1_FASTQC" ]; then
  UPLOAD_CMD="$UPLOAD_CMD --fastqcfile $R1_FASTQC"
fi
if [ -e "$R2_FASTQC" ]; then
  UPLOAD_CMD="$UPLOAD_CMD --fastqcfile $R2_FASTQC"
fi
if [ -e "$R3_FASTQC" ]; then
  UPLOAD_CMD="$UPLOAD_CMD --fastqcfile $R3_FASTQC"
fi
if [ -e "$R4_FASTQC" ]; then
  UPLOAD_CMD="$UPLOAD_CMD --fastqcfile $R4_FASTQC"
fi
eval $UPLOAD_CMD

#$APX bash $STAMPIPES/scripts/fastq/attachfiles.bash
# Inline contents of the script so I can add R3 and R4:
UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data_withr3.py
ATTACH_LANE="python3 $UPLOAD_SCRIPT --attach_file_contenttype SequencingData.flowcelllane --attach_file_object ${FLOWCELL_LANE_ID}"
$ATTACH_LANE --attach_directory ${FASTQ_DIR} --attach_file_purpose fastq-directory
if [ -e "$R1_FASTQC" ]; then
  $ATTACH_LANE --attach_file ${R1_FASTQC} --attach_file_type zip --attach_file_purpose fastqc-results-zip
fi
if [ -e "$R2_FASTQC" ]; then
  $ATTACH_LANE --attach_file ${R2_FASTQC} --attach_file_type zip --attach_file_purpose fastqc-results-zip
fi
if [ -e "$R3_FASTQC" ]; then
  $ATTACH_LANE --attach_file ${R3_FASTQC} --attach_file_type zip --attach_file_purpose fastqc-results-zip
fi
if [ -e "$R4_FASTQC" ]; then
  $ATTACH_LANE --attach_file ${R4_FASTQC} --attach_file_type zip --attach_file_purpose fastqc-results-zip
fi

echo "FINISH: "
date
