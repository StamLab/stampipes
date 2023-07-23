VERSION=1.2
OUT_DIR=output_$VERSION
mkdir -p "$OUT_DIR"

source "$MODULELOAD"
module purge
module load jdk nextflow openssl-dev

export REFDIR="$(dirname $GENOME_INDEX)"
export STARrefDir="$REFDIR/${STAR_DIR}"
export TARGET_BAM=Aligned.toTranscriptome.out.cram
export GENOME_BAM=Aligned.toGenome.out.cram
export NODUPS_BAM=Aligned.toGenome.noDups.cram
export TRIMS_R1=trims.R1.fastq.gz
export TRIMS_R2=trims.R2.fastq.gz

# Check for special UMTs
if [[ "$LIBRARY_KIT" == "SMARTer Stranded Total RNA-Seq Kit v3-Pico" ]] ; then
  UMI=True
  UMI_METHOD=takara-umt
fi

# record version
# TODO: Make more clear
cp "$STAMPIPES/version.json" "$OUT_DIR"

GENOME_BAM_FILES=$(sed 's/toTranscriptome/sortedByCoord/g' <<< "$BAM_FILES")

if [[ -e "$STARrefDir/all.genome.fa.gz" ]] ; then
  FASTAREF=$STARrefDir/all.genome.fa.gz
else
  FASTAREF=${GENOME_INDEX}.fa
fi

# Create params file
cat > agg_params.yaml <<PARAMS
id: "AG$AGGREGATION_ID"
umi: $UMI
annotation: "$ANNOTATION"
sequinsisomix: "$SEQUINS_ISO_MIX"
starrefdir: "$STARrefDir"
sequinsref: "$SEQUINS_REF"
kallistoindex: "$KALLISTO_INDEX"
neatmixa: "$NEAT_MIX_A"
flatref: "$FLAT_REF"
fastaref: "$FASTAREF"
outdir: "$OUT_DIR"

PARAMS

# Write bam files into parameter file
# Syntax is sensitive here, be careful when modifying.
echo "genomebams:" >> agg_params.yaml
for genomebam in $GENOME_BAM_FILES ; do
  echo "  - $genomebam" >> agg_params.yaml
done

echo "transcriptomebams:" >> agg_params.yaml
for transcriptomebam in $BAM_FILES ; do
  echo "  - $transcriptomebam" >> agg_params.yaml
done

set -e

# Run the code
NXF_VER=21.04.1 nextflow run \
  "$STAMPIPES/processes/rna-star/aggregation/cufflinks_featurecounts.nf" \
  -params-file agg_params.yaml \
  -profile cluster \
  -resume \
  -with-trace

# Upload results
(
  cd "$OUT_DIR"
  bash "$STAMPIPES/scripts/rna-star/aggregate/checkcomplete.sh"
  bash "$STAMPIPES/scripts/rna-star/aggregate/concat_metrics.sh"
  bash "$STAMPIPES/scripts/rna-star/aggregate/upload_counts.bash"
  bash "$STAMPIPES/scripts/rna-star/aggregate/attachfiles.sh"
)
