# Requires LIBRARY_NAME, AGGREGATION_ID and GENOME to be in the environment
# Removes any of the listed files that exist

echo "RESETTING AGGREGATION ${AGGREGATION_ID} FOR ${LIBRARY_NAME}"

files=( \
"adapter_counts.info" \
"Aligned.toGenome.out.bam" \
"Aligned.toGenome.out.bam.bai" \
"Aligned.toTranscriptome.out.bam" \
"chrNL.txt" \
"feature_counts.info" \
"feature_counts.txt" \
"feature_counts.txt.summary" \
"genes.fpkm_tracking" \
"isoforms.fpkm_tracking" \
"kallisto.log" \
"kallisto_adv.log" \ 
"metrics.info" \
"picard.CollectInsertSizes.txt" \
"picard.MarkDuplicates.txt" \
"picard.RnaSeqMetrics.txt" \
"ribosomal_counts.info" \
"rna_stats_summary.info" \
"Signal.Unique.both.starch.bgz" \
"Signal.Unique.both.starch.bgz.tbi" \
"Signal.UniqueMultiple.str-.bw" \
"Signal.UniqueMultiple.str+.bw" \
"Signal.UniqueMultiple.both.bw" \
"Signal.UniqueMultiple.str-.starch" \
"Signal.UniqueMultiple.str+.starch" \
"Signal.UniqueMultiple.both.starch" \
"Signal.Unique.str-.bw" \
"Signal.Unique.str+.bw" \
"Signal.Unique.both.bw" \
"Signal.Unique.str-.starch" \
"Signal.Unique.str+.starch" \
"Signal.Unique.both.starch" \
"skipped.gtf" \
"tagcounts.txt" \
"transcripts.gtf" \
"trims.R1.fastq.gz" \
"trims.R2.fastq.gz" \
)

dirs=( \
    "kallisto_output"       \
    "kallisto_output_adv"   \
)

for FILE in "${files[@]}"; do
    if [ -e "$FILE" ]; then
        echo "Removing $FILE"
        rm $FILE
    fi
done

for DIR in "${files[@]}"; do
    if [ -d $DIR ]; then
	echo "Removing Directory $DIR"
        rm -r $DIR
    fi
done

python3 $STAMPIPES/scripts/lims/upload_data.py --clear_aggregation_stats --aggregation_id ${AGGREGATION_ID}
