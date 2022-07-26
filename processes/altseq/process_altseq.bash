#!/bin/bash
# This script is copied by setup.sh to /net/seq/data2/flowcells/the_flowcell/

for var in FLOWCELL SEQUENCER_MOUNT ; do
  if [[ -z "$var" ]] ; then
    echo "Set env var $var"
    exit 2
  fi
done

set -eo pipefail

# TODO: Bump this before running for real.
version=0.9.1

cd "$(dirname "$0")"

outdir="output_$version"
sentinel_file="$outdir/process_complete.txt"

if [[ -e "$sentinel_file" && -z "$REDO_ALIGNMENT" ]] ; then
  echo "Processing already completed, exiting."
  echo "To force re-run, set the env var 'REDO_ALIGNMENT=True' or remove $sentinel_file"
  exit 0
fi

# Dependencies
source "$MODULELOAD"
module purge
module load jdk
module load nextflow
module load python/3.5.1

source "$PYTHON3_ACTIVATE"
source "$STAMPIPES/scripts/sentry/sentry-lib.bash"

# TODO: REDO_ALIGNMENT handling

# Set up sample config
sample_config=sample_config.tsv
python "$STAMPIPES"/scripts/lims/get_processing.py -f "$FLOWCELL"
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
  -profile docker \
  -resume \
  --input_directory "$SEQ_DIR"  \
  --sample_config_tsv "$sample_config" \
  --genome_dir "$GENOME_DIR" \
  --genome_fa "$GENOME_FA" \
  --barcode_whitelist "$BARCODE_WHITELIST" \
  --outdir "$outdir" \
  -ansi-log false


# Upload fastq metadata
python "$STAMPIPES/scripts/altseq/upload_fastq.py" \
  "$sample_config" \
  processing.json \
  --output_file_directory "$outdir"

if [[ ! -e "$sentinel_file" ]] ; then
  echo "{ completed_on: $(date -Iseconds) }" > "$sentinel_file"
fi
