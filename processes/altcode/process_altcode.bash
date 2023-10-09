#!/bin/bash
set -eo pipefail

version=1.0.0-alpha2
cd "$(dirname "$0")"

# Temporarily hardcoded!
R1_BARCODE_POS=78
R2_BARCODE_POS=48
R3_BARCODE_POS=10
R1_BARCODE_LEN=8
R2_BARCODE_LEN=8
R3_BARCODE_LEN=8

outdir="output_$version"
status_file="$outdir/status.json"

# TODO: improve REDO_ALIGNMENT handling - should we be manually removing the work dir?
if [[ -e "$status_file" && -z "$REDO_ALIGNMENT" ]] ; then
  # Check to see if the alignment is complete
  if jq -e '.completed_on' "$status_file" ; then
    echo "Processing already completed, exiting."
    echo "To force re-run, set the env var 'REDO_ALIGNMENT=True' or remove $status_file"
    exit 0
  fi
fi

# Dependencies
source "$MODULELOAD"
module purge
module load jdk/11.0.16
module load nextflow/22.04.3
module load python/3.5.1
module load apptainer/1.1.2

export NXF_VER=23.04.2

source "$PYTHON3_ACTIVATE"
source "$STAMPIPES/scripts/sentry/sentry-lib.bash"

WORKROOT=${WORKROOT:-/net/seq/scratch}
if ! [[ -d "$WORKROOT" ]] ; then
  echo "WORKROOT '$WORKROOT' does not exist, using '$PWD'"
  WORKROOT=$PWD
fi
WORKDIR=$WORKROOT/$USER/altseq/FC$FLOWCELL/work/


# Write parameter file
params=params.yaml
cat >$params <<PARAMS_YAML
outdir: "$outdir"
metadata: "pool_info.json"
r1: "$R1_FASTQ"
r2: "$R2_FASTQ"
genome_dir: "$STAR_GENOME_DIR"
genome_fasta: "$GENOME_FA"
barcode_r1_list: "$STAMPIPES_DATA/altcode/barcodes_r1_v2_96.txt"
barcode_r2_list: "$STAMPIPES_DATA/altcode/barcodes_r2.txt"
barcode_r3_list: "$STAMPIPES_DATA/altcode/barcodes_r3.txt"
barcode_r1_pos: $R1_BARCODE_POS
barcode_r2_pos: $R2_BARCODE_POS
barcode_r3_pos: $R3_BARCODE_POS
barcode_r1_len: $R1_BARCODE_LEN
barcode_r2_len: $R2_BARCODE_LEN
barcode_r3_len: $R3_BARCODE_LEN
umi_pos: 0
umi_len: 10
PARAMS_YAML

# Run :)
nextflow run \
    "$STAMPIPES/processes/altcode/altcode.nf" \
    -params-file "$params" \
    -ansi-log false \
    -with-trace \
    -work-dir "$WORKDIR" \
    -profile cluster \
    -resume

## Upload fastq metadata
#python "$STAMPIPES/scripts/altseq/upload_data.py" \
  #"$sample_config" \
  #processing.json \
  #--output_file_directory "$outdir"

# Create sentinel/status file
if [[ -e "$status_file" ]] ; then
  old_date=$(jq .completed_on <<< "$status_file")
  old_status_file=${status_file/json/$old_date}.json
  mv "$status_file" "$old_status_file"
fi

# TODO: What else do we want to capture here? It would be nice to at least
# capture the command used and relevant env vars
echo | jq . > "$status_file" <<EOF
  {
    "completed_on": "$(date -Iseconds)",
    "version": "$version"
  } 
EOF

