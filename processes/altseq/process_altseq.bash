#!/bin/bash
# This script is copied by setup.sh to /net/seq/data2/flowcells/the_flowcell/

for var in FLOWCELL SEQUENCER_MOUNT ; do
  if [[ -z "$var" ]] ; then
    echo "Set env var $var"
    exit 2
  fi
done

set -eo pipefail

version=1.0.0

cd "$(dirname "$0")"

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
module load jdk
module load nextflow
module load python/3.5.1

source "$PYTHON3_ACTIVATE"
source "$STAMPIPES/scripts/sentry/sentry-lib.bash"

# Set up sample config
sample_config=sample_config.tsv
# TODO: Re-enable after production time-out is fixed
# python "$STAMPIPES"/scripts/lims/get_processing.py -f "$FLOWCELL"
python "$STAMPIPES"/scripts/lims/create_altseq_sample_config.py processing.json --output "$sample_config"



SEQ_DIR=$(ls -d -1 ${SEQUENCER_MOUNT}/*$FLOWCELL* | head -n1)

GENOME_DIR=/net/seq/data2/projects/prime_seq/cell_ranger_ref/star_2.7.10_genome_2022_gencode.v39/
GENOME_FA=/net/seq/data2/projects/prime_seq/cell_ranger_ref/GRCh38-2022-Altius-gencode.v39-build/Homo_sapiens.GRCh38.dna.primary_assembly.fa.modified
BARCODE_WHITELIST=/net/seq/data2/projects/prime_seq/barcodes-combined.txt

# Run the pipeline
NXF_VER=21.10.6 nextflow \
  -c $STAMPIPES/nextflow.config \
  run "$STAMPIPES"/processes/altseq/altseq.nf \
  -with-trace \
  -ansi-log false \
  -profile docker,cluster \
  -resume \
  --input_directory "$SEQ_DIR"  \
  --sample_config_tsv "$sample_config" \
  --genome_dir "$GENOME_DIR" \
  --genome_fa "$GENOME_FA" \
  --barcode_whitelist "$BARCODE_WHITELIST" \
  --outdir "$outdir" \
  --skip_alignment


# Upload fastq metadata
python "$STAMPIPES/scripts/altseq/upload_data.py" \
  "$sample_config" \
  processing.json \
  --output_file_directory "$outdir"

# Create sentinel/status file
if [[ -e "$status_file" ]] ; then
  old_date=$(jq .completed_on << "$status_file")
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

# TODO: Enable this once we are running alignment
# nextflow clean -f .
