# Requires LIBRARY_NAME, AGGREGATION_ID, HOTSPOT2_DIR and GENOME to be in the environment
# Removes any of the listed files that exist

echo "RESETTING AGGREGATION ${AGGREGATION_ID} FOR ${LIBRARY_NAME}"

files=( \
    "${LIBRARY_NAME}.${GENOME}.sorted.bam" \ 
    "${LIBRARY_NAME}.${GENOME}.sorted.bam.bai" \ 
    "${LIBRARY_NAME}.tagcounts.txt" \
    "${LIBRARY_NAME}.${GENOME}.uniques.sorted.bam" \
    "${LIBRARY_NAME}.${GENOME}.uniques.sorted.bam.bai" \
    "${LIBRARY_NAME}.CollectInsertSizeMetrics.picard" \
    "${LIBRARY_NAME}.CollectInsertSizeMetrics.picard.pdf" \
    "${LIBRARY_NAME}.CollectInsertSizeMetrics.picard.info" \
    "${LIBRARY_NAME}.MarkDuplicates.picard" \
    "${LIBRARY_NAME}.75_20.${GENOME}.uniques-density.bed.starch" \
    "${LIBRARY_NAME}.75_20.${GENOME}.uniques-density.bed.starch.bgz" \
    "${LIBRARY_NAME}.75_20.${GENOME}.uniques-density.bed.starch.bgz.tbi" \
    "${LIBRARY_NAME}.75_20.${GENOME}.bw" \
    "${LIBRARY_NAME}.75_20.normalized.${GENOME}.bw" \
    "${LIBRARY_NAME}.75_20.normalized.${GENOME}.uniques-density.bed.starch" \ 
    "${LIBRARY_NAME}.75_20.normalized.${GENOME}.uniques-density.bed.starch.bgz" \
    "${LIBRARY_NAME}.75_20.normalized.${GENOME}.uniques-density.bed.starch.bgz.tbi" \
    "${LIBRARY_NAME}.${GENOME}.cuts.sorted.bed.starch" \
    "${LIBRARY_NAME}.${GENOME}.cutcounts.sorted.bed.starch" \
    "${LIBRARY_NAME}.${GENOME}.cutcounts.sorted.bed.starch.bgz" \
    "${LIBRARY_NAME}.${GENOME}.cutcounts.sorted.bed.starch.bgz.tbi" \
    "${LIBRARY_NAME}.${GENOME}.fragments.sorted.bed.starch" \
    "${LIBRARY_NAME}.${GENOME}.cutcounts.$READ_LENGTH.bw" \
    "${LIBRARY_NAME}.${GENOME}.R1.rand.uniques.sorted.spotdups.txt" \
    "${LIBRARY_NAME}.${GENOME}.R1.rand.uniques.sorted.spot.info" \
    "${LIBRARY_NAME}.${GENOME}.R1.rand.uniques.sorted.spot.out" \
    "${LIBRARY_NAME}.${GENOME}.uniques.sorted.hotspot2.info" \ 
    "${LIBRARY_NAME}.adaptercounts.txt" \
    "${LIBRARY_NAME}.${GENOME}.uniques.sorted.proxdist.info" \
    "${LIBRARY_NAME}.${GENOME}.rand.uniques.sorted.spotdups.txt" \
    "${LIBRARY_NAME}.${GENOME}.rand.uniques.sorted.spot.out" \
    "${LIBRARY_NAME}.${GENOME}.rand.uniques.sorted.spot.info" \
    "${LIBRARY_NAME}.uniques.duphist.txt" \
    "${LIBRARY_NAME}.uniques.preseq.txt" \
    "${LIBRARY_NAME}.uniques.preseq.targets.txt" \
)

for FILE in "${files[@]}"; do
    if [ -e "$FILE" ]; then
        echo "Removing $FILE"
        rm $FILE
    fi
done

if [[ -d "$HOTSPOT2_DIR" && "$HOTSPOT2_DIR" = "peaks_v2_1_1" ]]; then
    rm -r $HOTSPOT2_DIR
fi

python3 $STAMPIPES/scripts/lims/upload_data.py --clear_aggregation_stats --aggregation_id ${AGGREGATION_ID}
