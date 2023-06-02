# Uploads important files as attachments to an aggregation object on the LIMS
# Required environment names:
#  * STAMPIPES
#  * PEAK_CALLER
#  * AGGREGATION_ID

set -e -u -o pipefail

UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py

# All peaks files start with this:
PEAKS_PREFIX="peaks/nuclear"

ATTACH_AGGREGATION="python3 $UPLOAD_SCRIPT --attach_file_contenttype AggregationData.aggregation --attach_file_objectid ${AGGREGATION_ID}"

function attach_file() {
  local file=$1
  local purpose=$2
  local filetype=$3

  $ATTACH_AGGREGATION --attach_file "$file" \
    --attach_file_purpose "$purpose" \
    --attach_file_type "$filetype"
}

$ATTACH_AGGREGATION --attach_directory "$(readlink -f ..)" --attach_file_purpose aggregation-directory

attach_file  marked.cram       all-alignments-bam  cram
attach_file  marked.cram.crai  bam-index           cram-index

# Densities
attach_file  density.starch                density-bed-starch-windowed            starch
attach_file  density.bw                    density-bigwig-windowed                bigwig
attach_file  density.bgz                   density-tabix-bgz                      bgz
attach_file  normalized.density.starch     normalized-density-bed-starch          starch
attach_file  normalized.density.bw         normalized-density-bigwig-windowed     bigwig
attach_file  normalized.density.bgz        normalized-density-tabix-bgz           bgz
attach_file  mm_density.starch             mm-density-bed-starch-windowed         starch
attach_file  mm_density.bw                 mm-density-bigwig-windowed             bigwig
attach_file  normalized.mm_density.starch  normalized-mm-density-bed-starch       starch
attach_file  normalized.mm_density.bw      normalized-mm-density-bigwig-windowed  bigwig

# Cut counts
attach_file  cutcounts.bw       cutcounts-bw         bigwig
attach_file  cutcounts.starch   cutcounts-starch     starch
attach_file  cutcounts.bed.bgz  cutcounts-tabix-bgz  bgz

if [[ "$PEAK_CALLER" == hotspot2 ]] ; then
  # hotspot2 output
  attach_file  $PEAKS_PREFIX.allcalls.starch           hotspot-per-base                starch
  attach_file  $PEAKS_PREFIX.hotspots.fdr0.05.starch   hotspot-calls                   starch
  attach_file  $PEAKS_PREFIX.hotspots.fdr0.01.starch   hotspot-calls-1per              starch
  attach_file  $PEAKS_PREFIX.hotspots.fdr0.001.starch  hotspot-calls-point1per         starch
  attach_file  $PEAKS_PREFIX.peaks.fdr0.05.starch      hotspot-peaks                   starch
  attach_file  $PEAKS_PREFIX.peaks.fdr0.01.starch      hotspot-peaks-1per              starch
  attach_file  $PEAKS_PREFIX.peaks.fdr0.001.starch     hotspot-peaks-point1per         starch
  attach_file  $PEAKS_PREFIX.peaks.fdr0.001.bb         hotspot-peaks-point1per-bigbed  bigbed
fi

#TODO: We're effectively generating these twice, simplify
#$ATTACH_AGGREGATION --attach_file $PEAKS_PREFIX.cutcounts.starch --attach_file_purpose cutcounts-starch --attach_file_type starch
#$ATTACH_AGGREGATION --attach_file $PEAKS_PREFIX.density.starch --attach_file_purpose density-bed-starch-windowed --attach_file_type starch

python3 "$UPLOAD_SCRIPT" \
  --aggregation_id "${AGGREGATION_ID}" \
  --insertsfile CollectInsertSizeMetrics.picard

if [ -e MarkDuplicates.picard ]; then
  python3 "$UPLOAD_SCRIPT" \
    --aggregation_id "${AGGREGATION_ID}" \
    --dupsfile MarkDuplicates.picard
fi

if [ -e fragments.starch ]; then
  attach_file fragments.starch fragments-starch starch
fi
