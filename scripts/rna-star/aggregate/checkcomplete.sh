# Requires SAMPLE_NAME and GENOME to be in the environment
# Checks that important files exist and are not size 0

EXIT=0

files=( \
    "Aligned.toGenome.out.bam" \
    "Aligned.toTranscriptome.out.bam" \
    "feature_counts.txt" \
    "genes.fpkm_tracking" \
    "isoforms.fpkm_tracking" \
    "Signal.UniqueMultiple.str-.bw" \
    "Signal.UniqueMultiple.str+.bw" \
)

for FILE in "${files[@]}"; do
    if [ ! -s $FILE ]; then
	echo "Missing $FILE"
	EXIT=1
    fi
done

if [[ $EXIT -ne 1 ]]; then
    python3 "$STAMPIPES/scripts/lims/upload_data.py" --aggregation_id ${AGGREGATION_ID} --complete_aggregation
fi

exit $EXIT
