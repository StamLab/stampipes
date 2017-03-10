# Requires SAMPLE_NAME and GENOME to be in the environment
# Checks that important files exist and are not size 0

EXIT=0

files=( \
    "${SAMPLE_NAME}.adapters.txt" \
    "${SAMPLE_NAME}.versions.txt" \
    "Aligned.sortedByCoord.out.bam" \
    "Aligned.toTranscriptome.out.bam" \
)
    #"Quant.genes.results" \
    #"Quant.isoforms.results" \
    #"Quant.pdf" \
    #"Signal.Unique.strand+.bw" \
    #"Signal.Unique.strand-.bw" \
    #"Signal.UniqueMultiple.strand+.bw" \
    #"Signal.UniqueMultiple.strand-.bw" \

for FILE in "${files[@]}"; do
if [ ! -s $FILE ]; then
    echo "Missing $FILE"
    EXIT=1
fi
done

exit $EXIT
