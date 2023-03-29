# Uploads important files as attachments to an aggregation object on the LIMS
# Required environment names:
#  * STAMPIPES
#  * AGGREGATION_ID

UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py

function attach_aggregation () {
    python3 "$UPLOAD_SCRIPT" \
        --attach_file_contenttype AggregationData.aggregation \
        --attach_file_objectid "$AGGREGATION_ID" \
        "$@"
}
function attach_agg_file () {
    local filename=$1
    local purpose=$2
    local filetype=$3
    attach_aggregation \
        --attach_file "$filename"  \
        --attach_file_purpose "$purpose" \
        --attach_file_type "$filetype"
}


attach_aggregation --attach_directory "$PWD/.." --attach_file_purpose aggregation-directory

# alignments
attach_agg_file dupsmarked.cram all-alignments-bam cram
attach_agg_file merged.transcriptome.cram transcriptome-alignments cram

# cufflinks
attach_agg_file genes.fpkm_tracking gene-expression-levels plaintext
attach_agg_file isoforms.fpkm_tracking isoform-expression-levels plaintext

# kallisto
attach_agg_file kallisto_output/abundance.tsv kallisto-abundance tabtext

# featureCounts
attach_agg_file feature_counts.txt featureCounts-output plaintext

# densities
attach_agg_file Signal.UniqueMultiple.str+.bw pos-coverage-bigwig bigwig
attach_agg_file Signal.UniqueMultiple.str-.bw neg-coverage-bigwig bigwig
attach_agg_file Signal.UniqueMultiple.both.bw all-coverage-bigwig bigwig

# trimmed fastqs
# attach_agg_file trims.R1.fastq.gz r1-fastq-trimmed gzipped-fastq
# attach_agg_file trims.R2.fastq.gz r2-fastq-trimmed gzipped-fastq

# picard uploads
python3 "$UPLOAD_SCRIPT" --aggregation_id "$AGGREGATION_ID" --insertsfile picard.CollectInsertSizes.txt
python3 "$UPLOAD_SCRIPT" --aggregation_id "$AGGREGATION_ID" --dupsfile picard.MarkDuplicates.txt
