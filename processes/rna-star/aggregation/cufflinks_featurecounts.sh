source $MODULELOAD
module load samtools/1.7
module load gcc/4.7.2     # R dependent on this
module load R/3.2.5
module load STAR/2.4.2a   # Just for densities
module load bedops/2.4.19
module load subread/1.5.1 # for featureCounts
module load cufflinks/2.2.1 # for cuffLinks
module load anaquin/2.0.1
module load kallisto/0.43.1
module load htslib/1.6.0
module load jdk/1.8.0_92
module load picard/2.8.1
source "$PYTHON3_ACTIVATE"
module load python/2.7.11
module load bowtie/1.0.0
module load stringtie/1.3.4d

export REFDIR="$(dirname $GENOME_INDEX)"
export STARrefDir="$REFDIR/${STAR_DIR}"
export TARGET_BAM=Aligned.toTranscriptome.out.bam
export GENOME_BAM=Aligned.toGenome.out.bam
export NODUPS_BAM=Aligned.toGenome.noDups.bam
export TRIMS_R1=trims.R1.fastq.gz
export TRIMS_R2=trims.R2.fastq.gz

if [[ -n "$UMI_METHOD" ]] ; then
  export HAVE_UMI=1
  export BAM_TO_USE=$NODUPS_BAM
else 
  export HAVE_UMI=0
  export BAM_TO_USE=$GENOME_BAM
fi

# This will be set to the SLURM job IDs to wait for if we are removing dups
WAIT_FOR_DUPS=

function make_dependency_param() {
  # input: List of jobids, like  `,1234,1235,1236`
  # output: Suitable for SLURM input, like `--dependency=afterok:1234:1235:1236`
  local jobids=$1
  echo "$jobids" | sed -e 's/^,/--dependency=afterany:/;s/,/:/g'
}

if [ -n "$REDO_AGGREGATION" ]; then
    bash $STAMPIPES/scripts/rna-star/aggregate/reset.bash
fi

# record version
cp $STAMPIPES/version.json .

# create proper merged BAM
numbam=$(wc -w <<< $BAM_FILES)
# Temporary
if [ ! -s "$TARGET_BAM" ] ; then
  if [ $numbam -eq 1 ] ; then
    cp "$BAM_FILES" "$TARGET_BAM"
  else
    samtools merge -n -f "$TARGET_BAM" $BAM_FILES
  fi
fi
if [ ! -s "$GENOME_BAM" ] ; then
  GENOME_BAM_FILES=$(sed 's/toTranscriptome/sortedByCoord/g' <<< "$BAM_FILES")
  $STAMPIPES/scripts/tophat/merge_or_copy_bam.sh "$GENOME_BAM" $GENOME_BAM_FILES
  samtools index "$GENOME_BAM"
fi

# create merged fastqs
if [ ! -s "$TRIMS_R1" ] ; then
  samtools sort -n "$GENOME_BAM" -o sorted.bam
  samtools fastq --threads 2 sorted.bam -1 "$TRIMS_R1" -2 "$TRIMS_R2"
  rm sorted.bam
fi

density_job=.AG${AGGREGATION_ID}.star_den
cufflinks_job=.AG${AGGREGATION_ID}.star_cuff
kallisto_job=.AG${AGGREGATION_ID}.kallisto
kallisto_adv_job=.AG${AGGREGATION_ID}.kallisto_adv
complete_job=.AG${AGGREGATION_ID}.complete
fcounts_job=.AG${AGGREGATION_ID}.star_fcounts
stringtie_job=.AG${AGGREGATION_ID}.stringtie
picard_job=.AG${AGGREGATION_ID}.picard
dupes_job=.AG${AGGREGATION_ID}.dupes
adaptercounts_job=.AG${AGGREGATION_ID}.adapter
anaquin_job=.AG${AGGREGATION_ID}.anaquinsub
bow_rRNA_job=.AG${AGGREGATION_ID}.rRNA

# density information, convoluted, can clean up so we skip a lot of these steps
if [ ! -s "Signal.Unique.both.starch.bgz.tbi" ] ; then
    jobid=$(sbatch --export=ALL -J "$density_job" -o "$density_job.o%A" -e "$density_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__DEN__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START DENSITY: "
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

    # Write starch and bigwig to .tmp files
    function convertBedGraph(){
      in="\$1"
      base="\$2"
      chrom="\$in.onlyChr.bg"
      grep '^chr' "\$in" | sort -k1,1 -k2,2n > \$chrom
      bedGraphToBigWig "\$chrom" chrNL.txt "\$base.bw.tmp"
      starch "\$chrom" > "\$base.starch.tmp"
    }

    mkdir -p \$TMPDIR/Signal

    echo STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix \$TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir \$TMPDIR/STAR
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Unstranded --outFileNamePrefix \$TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir \$TMPDIR/STAR
    mv \$TMPDIR/Signal/Signal.UniqueMultiple.str1.out.bg \$TMPDIR/Signal/Signal.UniqueMultiple.unstranded.out.bg
    mv \$TMPDIR/Signal/Signal.Unique.str1.out.bg \$TMPDIR/Signal/Signal.Unique.unstranded.out.bg
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix \$TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir \$TMPDIR/STAR

    grep '^chr' $STARrefDir/chrNameLength.txt > chrNL.txt

    convertBedGraph \$TMPDIR/Signal/Signal.Unique.str1.out.bg         Signal.Unique.str-
    convertBedGraph \$TMPDIR/Signal/Signal.Unique.str2.out.bg         Signal.Unique.str+
    convertBedGraph \$TMPDIR/Signal/Signal.Unique.unstranded.out.bg         Signal.Unique.both
    convertBedGraph \$TMPDIR/Signal/Signal.UniqueMultiple.str1.out.bg Signal.UniqueMultiple.str-
    convertBedGraph \$TMPDIR/Signal/Signal.UniqueMultiple.str2.out.bg Signal.UniqueMultiple.str+
    convertBedGraph \$TMPDIR/Signal/Signal.UniqueMultiple.unstranded.out.bg Signal.UniqueMultiple.both

    for i in Signal*.tmp ; do
      mv \$i \${i/.tmp/}
    done

    unstarch Signal.Unique.both.starch | bgzip > Signal.Unique.both.starch.bgz
    tabix -p bed Signal.Unique.both.starch.bgz

rm -rf "\$TMPDIR"

echo "FINISH DENSITY: "
date

__DEN__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# cufflinks
if [ ! -s "isoforms.fpkm_tracking" ] ; then
    jobid=$(sbatch --export=ALL -J "$cufflinks_job" -o "$cufflinks_job.o%A" -e "$cufflinks_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__CUFF__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START CUFFLINKS: "
date

    CUFF=cufflinks
    CUFF_COMMON="--no-update-check --library-type fr-firststrand"

    \$CUFF \$CUFF_COMMON --GTF $ANNOTATION $GENOME_BAM

    # delete duplicate rows and sort into identical orders across all samples
    Rscript $STAMPIPES/scripts/rna-star/aggregate/dedupe_sort_cuffout.Rscript genes.fpkm_tracking
    Rscript $STAMPIPES/scripts/rna-star/aggregate/dedupe_sort_cuffout.Rscript isoforms.fpkm_tracking
    mv genes.fpkm_tracking.sort genes.fpkm_tracking
    mv isoforms.fpkm_tracking.sort isoforms.fpkm_tracking

    # quantification with anaquin Rna Expression
    if [ -n "$SEQUINS_REF" ] ; then
        anaquin RnaExpression -o anaquin_cufflinks -rmix $SEQUINS_ISO_MIX -usequin transcripts.gtf -mix A || (echo "NA" > anaquin_cufflinks/RnaExpression_genes.tsv && echo "NA" > anaquin_cufflinks/RnaExpression_isoforms.tsv && echo "NA" > anaquin_cufflinks/RnaExpression_summary.stats)
    fi

echo "FINISH CUFFLINKS: "
date

__CUFF__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# stringtie
if [ ! -s "stringtie_rf/abundances.txt" ] ; then
    jobid=$(sbatch --export=ALL -J "$stringtie_job" -o "$stringtie_job.o%A" -e "$stringtie_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__ST__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START STRINGTIE: "
date

stringtie $GENOME_BAM --rf -p 4 -G $ANNOTATION -o stringtie_rf/transcripts.gtf -A stringtie_rf/abundances.txt

echo "FINISH STRINGTIE: "
date

__ST__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# featureCounts
if [ ! -s "feature_counts.txt" ] ; then
    jobid=$(sbatch --export=ALL -J "$fcounts_job" -o "$fcounts_job.o%A" -e "$fcounts_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__FC__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START FEATURECOUNTS: "
date

    FCOUNTS=featureCounts
    FCOUNTS_COMMON="--primary -B -C -p -P --fracOverlap .5 -s 2 -D 10000"
    \$FCOUNTS \$FCOUNTS_COMMON -t 'exon' -g 'gene_id' -a $ANNOTATION -o feature_counts.txt $GENOME_BAM

echo "FINISH FEATURECOUNTS: "
date

__FC__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# kallisto
if [ ! -s "kallisto_output/abundance.tsv" ] ; then
    jobid=$(sbatch --export=ALL -J "$kallisto_job" -o "$kallisto_job.o%A" -e "$kallisto_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__KALLISTO__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START KALLISTO: "
date

kallisto quant -i $KALLISTO_INDEX -o kallisto_output $TRIMS_R1 $TRIMS_R2 2> kallisto.log

if [ -n "\$SEQUINS_REF" ] ; then
    anaquin RnaExpression -o anaquin_kallisto -rmix $SEQUINS_ISO_MIX -usequin kallisto_output/abundance.tsv -mix A || (echo "NA" > anaquin_kallisto/RnaExpression_genes.tsv && echo "NA" > anaquin_kallisto/RnaExpression_isoforms.tsv && echo "NA" > anaquin_kallisto/RnaExpression_summary.stats)
fi

echo "FINISH KALLISTO: "
date

__KALLISTO__
)
    PROCESSING="$PROCESSING,$jobid"
fi

# kallisto advanced
if [ ! -s "kallisto_output_adv/abundance.tsv" ] ; then
    jobid=$(sbatch --export=ALL -J "$kallisto_adv_job" -o "$kallisto_adv_job.o%A" -e "$kallisto_adv_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=64000 --parsable --oversubscribe <<__KALLISTO__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START KALLISTO ADV: "
date

kallisto quant --bias -b 100 --rf-stranded -i $KALLISTO_INDEX -o kallisto_output_adv $TRIMS_R1 $TRIMS_R2 2> kallisto_adv.log

if [ -n "\$SEQUINS_REF" ] ; then
    anaquin RnaExpression -o anaquin_kallisto_adv -rmix $SEQUINS_ISO_MIX -usequin kallisto_output_adv/abundance.tsv -mix A || (echo "NA" > anaquin_kallisto_adv/RnaExpression_genes.tsv && echo "NA" > anaquin_kallisto_adv/RnaExpression_isoforms.tsv && echo "NA" > anaquin_kallisto_adv/RnaExpression_summary.stats)
fi

echo "FINISH KALLISTO ADV: "
date

__KALLISTO__
)
    PROCESSING="$PROCESSING,$jobid"
fi

# picard
if [[ ! -s "picard.CollectInsertSizes.txt" || ! -s "picard.RnaSeqMetrics.txt" || ! -s "rna_stats_summary.info" ]] ; then
    jobid=$(sbatch --export=ALL -J "$picard_job" -o "$picard_job.o%A" -e "$picard_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__PIC__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START PICARD: "
date

picard CollectInsertSizeMetrics INPUT=$GENOME_BAM OUTPUT=picard.CollectInsertSizes.txt HISTOGRAM_FILE=/dev/null
picard CollectRnaSeqMetrics INPUT=$GENOME_BAM OUTPUT=picard.RnaSeqMetrics.txt REF_FLAT=$FLAT_REF STRAND_SPECIFICITY="SECOND_READ_TRANSCRIPTION_STRAND"

cat picard.RnaSeqMetrics.txt | grep -A 1 "METRICS CLASS" | sed 1d | tr '\t' '\n' > rna_stats_summary.info
cat picard.RnaSeqMetrics.txt | grep -A 2 "METRICS CLASS" | sed 1d | sed 1d | tr '\t' '\n' | paste rna_stats_summary.info - > tmp.txt && mv tmp.txt rna_stats_summary.info

echo "FINISH PICARD: "
date

__PIC__
)
   PROCESSING="$PROCESSING,$jobid"
fi

# dups
if [ ! -s "picard.MarkDuplicates.txt" ] ; then
    jobid=$(sbatch --export=ALL -J "$dupes_job" -o "$dupes_job.o%A" -e "$dupes_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__TC__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START DUPLICATES: "
date

if [[ "$HAVE_UMI" == 1 ]] ; then
  # Add Mate Cigar information; required by UMI-aware MarkDuplicates
  picard RevertOriginalBaseQualitiesAndAddMateCigar \
    "INPUT=$GENOME_BAM OUTPUT=cigar.bam \
    VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

  # Remove non-primary reads (also required)
  # TODO: Do we need -f2 ?
  samtools view -f2 -F256 cigar.bam -o cigar_no_supp.bam
  picard UmiAwareMarkDuplicatesWithMateCigar \
    INPUT=cigar_no_supp.bam \
    OUTPUT=$NODUPS_BAM \
    METRICS_FILE=picard.MarkDuplicates.txt \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' \
    UMI_TAG_NAME=RX \
    UMI_METRICS=umi_metrics.txt \
    REMOVE_DUPLICATES=true

  # Remove intermediate files
  rm cigar.bam cigar_no_supp.bam
else
  # No UMI info, proceed as normal
  picard MarkDuplicates \
    INPUT=$GENOME_BAM \
    OUTPUT=/dev/null \
    METRICS_FILE=picard.MarkDuplicates.txt \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
fi


echo "FINISH DUPLICATES: "
date

__TC__
)
  PROCESSING="$PROCESSING,$jobid"
  if [[ $HAVE_UMI == 1 ]] ; then
    WAIT_FOR_DUPS=$(make_dependency_param ",$jobid")
  fi
fi


# adapter counts
if [ ! -s "adapter_counts.info" ] ; then
    jobid=$(sbatch --export=ALL -J "$adaptercounts_job" -o "$adaptercounts_job.o%A" -e "$adaptercounts_job.e%A" $WAIT_FOR_DUPS --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: adapter counts"
date

adaptercount=\$(bash "$STAMPIPES/scripts/bam/count_adapters.sh" "$BAM_TO_USE")
if [ -n \$adapter_count ]; then
        echo -e "adapter\t\$adaptercount" > "adapter_counts.info"
fi

echo "FINISH: adapter counts"
date

__SCRIPT__
)
    PROCESSING="$PROCESSING,$jobid"
fi

# rRNA
if [ ! -s "ribosomal_counts.info" ] ; then
    jobid=$(sbatch --export=ALL -J "$bow_rRNA_job" -o "$bow_rRNA_job.o%A" -e "$bow_rRNA_job.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

echo "Hostname: "
hostname

echo "START BOWTIE: "
date

zcat -f $TRIMS_R1 > \$TMPDIR/trims.R1.fastq
zcat -f $TRIMS_R2 > \$TMPDIR/trims.R2.fastq
num_reads=\$(bowtie -n 3 -e 140 /net/seq/data/genomes/human/GRCh38/noalts/contamination/hg_rRNA -1 \$TMPDIR/trims.R1.fastq -2 \$TMPDIR/trims.R2.fastq | wc -l )
if [ -n \$num_reads ]; then
    echo -e "ribosomal-RNA\t\$num_reads" > "ribosomal_counts.info"
fi

rm -rf "\$TMPDIR"

echo "FINISH BOWTIE: "
date

__SCRIPT__
)
    PROCESSING="$PROCESSING,$jobid"
fi

# this correctly fails out if there are no sequins alignments found
# this creates a weird fail state if there are sequins alignments found but not enough to do the subsampling, not sure how to fix this other than to change anaquin to fail as it were the above case
# downstream uploads of information still occur
if [[ ! -s "anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info" && -n "$SEQUINS_REF" ]] ; then
    jobid=$(sbatch --export=ALL -J "$anaquin_job" -o "$anaquin_job.o%A" -e "$anaquin_job.e%A" --partition=$QUEUE $WAIT_FOR_DUPS --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__ANAQUIN__
#!/bin/bash

set -x -e -o pipefail

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

echo "Hostname: "
hostname

echo "START ANAQUIN: "
date

# calculate anaquin align stats on full alignment
anaquin RnaAlign -rgtf $SEQUINS_REF -usequin $BAM_TO_USE -o anaquin_star
bash $STAMPIPES/scripts/rna-star/aggregate/anaquin_rnaalign_stats.bash anaquin_star/RnaAlign_summary.stats anaquin_star/RnaAlign_summary.stats.info

# create subsample
DILUTION="0.0001"
anaquin RnaSubsample -method \$DILUTION -usequin $BAM_TO_USE -o anaquin_subsample | samtools view -bS - > \$TMPDIR/temp_subsample.bam

# calculate anaquin align stats on subset alignment
anaquin RnaAlign -rgtf $SEQUINS_REF -usequin \$TMPDIR/temp_subsample.bam -o anaquin_subsample/anaquin_star
bash $STAMPIPES/scripts/rna-star/aggregate/anaquin_rnaalign_stats.bash anaquin_subsample/anaquin_star/RnaAlign_summary.stats anaquin_subsample/anaquin_star/RnaAlign_summary.stats.info

# turn subset alignment to fastq
samtools sort -n -o \$TMPDIR/temp_subsample.sorted.bam \$TMPDIR/temp_subsample.bam
samtools view -uf64 \$TMPDIR/temp_subsample.sorted.bam | samtools fastq - > \$TMPDIR/subsample.fq1
samtools view -uf128 \$TMPDIR/temp_subsample.sorted.bam | samtools fastq - > \$TMPDIR/subsample.fq2

# call kallisto on subsampled fastqs
kallisto quant -i $KALLISTO_INDEX -o anaquin_subsample/kallisto_output \$TMPDIR/subsample.fq1 \$TMPDIR/subsample.fq2

# calculate kallisto stats on subsampled alignment
anaquin RnaExpression -o anaquin_subsample/anaquin_kallisto -rmix $SEQUINS_ISO_MIX -usequin anaquin_subsample/kallisto_output/abundance.tsv -mix A || (echo "NA" > anaquin_subsample/anaquin_kallisto/RnaExpression_genes.tsv && echo "NA" > anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.tsv && echo "NA" > anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats)
bash $STAMPIPES/scripts/rna-star/aggregate/anaquin_rnaexp_stats.bash anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info
bash $STAMPIPES/scripts/rna-star/aggregate/anaquin_neat_comparison.bash anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.tsv anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv $NEAT_MIX_A

rm -rf "\$TMPDIR"

echo "FINISH ANAQUIN: "
date

__ANAQUIN__
)
    PROCESSING="$PROCESSING,$jobid"
fi

# complete
dependencies_full=$(echo $PROCESSING | sed -e 's/,/,afterany:/g' | sed -e 's/^,afterany/--dependency=afterok/g')
sbatch --export=ALL -J "$complete_job" -o "$complete_job.o%A" -e "$complete_job.e%A" $dependencies_full --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=8000 --parsable --oversubscribe <<__COMPLETE__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START COMPLETE: "
date

# manually delete extremely large anaquin star log files (cannot be suppressed)
if [[ -s "anaquin_star/anaquin.log" ]] ; then
   rm anaquin_star/anaquin.log
fi
if [[ -s "anaquin_subsample/anaquin_star/anaquin.log" ]] ; then
   rm anaquin_subsample/anaquin_star/anaquin.log
fi

bash $STAMPIPES/scripts/rna-star/aggregate/checkcomplete.sh
bash $STAMPIPES/scripts/rna-star/aggregate/concat_metrics.sh
bash $STAMPIPES/scripts/rna-star/aggregate/upload_counts.bash
bash $STAMPIPES/scripts/rna-star/aggregate/attachfiles.sh

echo "FINISH COMPLETE: "
date

__COMPLETE__
