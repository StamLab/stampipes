# Requires LIBRARY_NAME, AGGREGATION_ID and GENOME to be in the environment
# Removes any of the listed files that exist

echo "RESETTING AGGREGATION ${AGGREGATION_ID} FOR ${LIBRARY_NAME}"

files=( \
"Aligned.toGenome.out.bam" \
"Aligned.toGenome.out.bam.bai" \
"Aligned.toTranscriptome.out.bam" \
"feature_counts.txt" \
"feature_counts.txt.summary" \
"genes.fpkm_tracking" \
"isoforms.fpkm_tracking" \
"Signal.UniqueMultiple.str-.bw" \
"Signal.UniqueMultiple.str+.bw" \
"Signal.UniqueMultiple.str-.starch" \
"Signal.UniqueMultiple.str+.starch" \
"Signal.Unique.str-.bw" \
"Signal.Unique.str+.bw" \
"Signal.Unique.str-.starch" \
"Signal.Unique.str+.starch" \
"skipped.gtf" \
"transcripts.gtf" \
)

for FILE in "${files[@]}"; do
    if [ -e "$FILE" ]; then
        echo "Removing $FILE"
        rm $FILE
    fi
done

python3 $STAMPIPES/scripts/lims/upload_data.py --clear_aggregation_stats --aggregation_id ${AGGREGATION_ID}