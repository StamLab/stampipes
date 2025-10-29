#!/bin/bash
set -euo pipefail

# TODO: Move this logic to another script probably - running this to make sure python can run upload_data.py
CLUSTER_NAME=$(scontrol show config | awk '$1 == "ClusterName" {print $3}')
if [[ "$CLUSTER_NAME" == "altius-gene" ]]; then
  module load apptainer/1.3.3
  echo "# Using apptainer"
  export APX="apptainer exec --bind '/net/seq/data2/flowcells,$STAMPIPES' '$STAMPIPES/containers/fastq/fastq.sif'"
else
  echo "# Not using apptainer"
  export APX=
fi
# END TODO

# Constants (do not depend on sample metadata)
export NXF_VER=25.04.6
NFCORE_RNASEQ_VER="3.21.0"
PIPELINE_VER=1.0-alpha
THIS_DIR=$(readlink -f "$(dirname "$0")")   # Will resolve to the location of the running bash script - e.g: /path/to/Sample_DS12345/align.../
OUT_DIR="$THIS_DIR/output_$PIPELINE_VER"
PROCESS_CONFIG="$STAMPIPES/processes/nf-core/rnaseq/nextflow.config"
SAMPLESHEET="$THIS_DIR/samplesheet.csv"

# Variables based on metadata vars (written by alignprocess.py)
STRANDEDNESS="reverse"
# Example:
#REFDIR=$(dirname "$BWAINDEX")
#export STARdir="$REFDIR/$STAR_DIR"
#bash "$STAMPIPES/scripts/versions.bash" &> "$VERSION_FILE"

mkdir -p "$OUT_DIR"
cd "$THIS_DIR"

# TODO: Remove this example of matching on library kit
# Copied from existing RNA star pipeline
## Check for special UMTs
#if [[ "$LIBRARY_KIT" == "Smarter Stranded Total v3 Pico RNASeq with RNA Isolation" ]] ; then
#  UMI=True
#  UMI_METHOD=takara-umt
#fi
#
#ADAPTER_FLAG=--use_fastp
#if [[ "$LIBRARY_KIT" == "Agilent Sure Select XT HS2" ]] ; then
#  ADAPTER_FLAG=--use_agent
#fi

# If we have a module system; try to load nextflow modules
set +u
if [[ -n "$MODULELOAD" || -n "$MODULEPATH" ]] ; then
  module load nextflow/25.04.6
fi
set -u

# Track versions
VERSION_FILE="$OUT_DIR/$SAMPLE_NAME.versions.txt"

cat > "$VERSION_FILE" <<- EOM
pipeline: $PIPELINE_VER
nextflow: $NXF_VER
nf-core/rnaseq: $NFCORE_RNASEQ_VER
EOM



# Set up sample sheet
cat > "$SAMPLESHEET" <<- EOM
sample,fastq_1,fastq_2,strandedness
$SAMPLE_NAME,$FASTQ_R1,$FASTQ_R2,$STRANDEDNESS
EOM

# Let LIMS know the alignment is starting
$APX python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  -a "$LIMS_API_URL"             \
  -t "$LIMS_API_TOKEN"           \
  --alignment_id "$ALIGNMENT_ID" \
  --start_alignment_progress     \
  --version_file "$VERSION_FILE"

# Run the pipeline
nextflow run "nf-core/rnaseq" \
  -r "$NFCORE_RNASEQ_VER" \
  -c "$PROCESS_CONFIG" \
  -profile cluster \
  -resume \
  "$@"

# TODO: Add post-processing steps here
# Check for completeness and upload files.
# set -e
# bash "$STAMPIPES/scripts/rna-star/checkcomplete.bash"
# bash "$STAMPIPES/scripts/rna-star/attachfiles.sh"
# 
# # Signal all-complete
# $APX python3 "$STAMPIPES/scripts/lims/upload_data.py" \
#   -a "$LIMS_API_URL"             \
#   -t "$LIMS_API_TOKEN"           \
#   --alignment_id "$ALIGNMENT_ID" \
#   --finish_alignment
