#!/bin/bash

version=2.6.0
export NXF_VER=18.10.1  # The version of nextflow to run. 18.10.1 includes conda

cd "$(dirname "$0")"

source "$MODULELOAD"
module purge
module load python/3.5.1
module load anaconda/2.1.0-dev

source "$PYTHON3_ACTIVATE"
source "$STAMPIPES/scripts/sentry/sentry-lib.bash"

if [[ $(wc -w <<< "$BAM_FILES") -gt 1 ]] ; then
  bamfiles="$(sed 's/\s\+/,/g' <<< "$BAM_FILES")"
else
  bamfiles=$BAM_FILES
fi
outdir="output_${version}"
workdir="work"

WIN=75
BINI=20


export MAPPABLE_REGIONS=${MAPPABLE_REGIONS:-$GENOME_INDEX.K${READ_LENGTH}.mappable_only.bed}
export CHROM_SIZES=${CHROM_SIZES:-$GENOME_INDEX.chrom_sizes.bed}
export CENTER_SITES=${CENTER_SITES:-$GENOME_INDEX.K${READ_LENGTH}.center_sites.n100.nuclear.starch}
export NUCLEAR_CHR=${NUCLEAR_CHR:-$GENOME_INDEX.nuclear.txt}
export CHROM_BUCKET=$STAMPIPES_DATA/densities/chrom-buckets.$GENOME.${WIN}_${BINI}.bed.starch

HOTSPOT_INDEX=${HOTSPOT_INDEX:-.}

# Remove old stuff if necessary
if [[ -n "$REDO_AGGREGATION" ]] ; then
  rm -rf "$outdir"
  rm -rf "$workdir"
  echo \
  python3 "$STAMPIPES/scripts/lims/upload_data.py" \
    --clear_aggregation_stats \
    --aggregation_id "$AGGREGATION_ID"
fi

# Tell LIMS we're starting alignment
python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  -a "$LIMS_API_URL" \
  -t "$LIMS_API_TOKEN" \
  --aggregation_id "$AGGREGATION_ID" \
  --start_aggregation

# Set peak_caller to macs2 for Chip-seq data
if [[ -z "$PEAK_CALLER" ]] ; then
  if [[ "$ASSAY_CATEGORY" == chipseq || "$ASSAY_CATEGORY" == fiberseq ]] ; then
    PEAK_CALLER=macs2
  else
    PEAK_CALLER=hotspot2
  fi
fi

# Run the whole process
(
  module purge
  module load jdk nextflow
  nextflow run \
    "$STAMPIPES/processes/bwa/aggregate/basic.nf" \
    -c "$STAMPIPES/nextflow.config" \
    -w "$workdir" \
    --bams "$bamfiles" \
    --genome "$GENOME_INDEX" \
    --mappable "$MAPPABLE_REGIONS" \
    --chrom_sizes "$CHROM_SIZES" \
    --peakcaller "$PEAK_CALLER" \
    --centers "$CENTER_SITES" \
    --chrom_bucket "$CHROM_BUCKET" \
    --hotspot_index "$HOTSPOT_INDEX" \
    --hotspot_id "AG$AGGREGATION_ID" \
    --bias "$STAMPIPES_DATA/footprints/vierstra_et_al.txt" \
    --UMI="$UMI" \
    --outdir "$outdir" \
    --threads 3 \
    -profile cluster,modules \
    -with-report nextflow.report.html \
    -with-dag nextflow.flowchart.html \
    -with-timeline nextflow.timeline.html \
    -resume
)

export PEAK_CALLER
( cd "$outdir" \
    && bash "$STAMPIPES/scripts/bwa/aggregate/basic/attachfiles_nextflow.bash" \
    && bash "$STAMPIPES/scripts/bwa/aggregate/basic/uploadcounts_nextflow.bash"
)
# Upload results
python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  -a "$LIMS_API_URL" \
  -t "$LIMS_API_TOKEN" \
  --aggregation_id "$AGGREGATION_ID" \
  --complete_aggregation
