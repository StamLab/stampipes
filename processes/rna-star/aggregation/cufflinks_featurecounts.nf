nextflow.enable.dsl=2

params.outdir = "output"
params.genomebams = []
params.transcriptomebams = []
params.umi = false

params.annotation = null
params.sequinsisomix = null
params.starrefdir = null
params.kallistoindex = null
params.sequinsref = null
params.neatmixa = null
params.flatref = null


workflow {
  RNA_AGG()
}

workflow RNA_AGG {

  // Grab params
  transcriptome_bams = params.transcriptomebams.collect {file it}
  genome_bams = params.genomebams.collect {file it}
  annotation = file(params.annotation)
  sequins_iso_mix = file(params.sequinsisomix)
  star_ref_files = file("${params.starrefdir}/*")

  sequins_ref = file(params.sequinsref)
  kallisto_index = file(params.kallistoindex)
  neat_mix_A = file(params.neatmixa)
  flat_ref = file(params.flatref)


  // Prep work
  merge_transcriptome_bam(transcriptome_bams)
  merge_genome_bam(genome_bams)

  // Remove UMI if present
  if (params.umi) {
    remove_duplicate_reads(merge_genome_bam.out.bam)
    bam_to_use = remove_duplicate_reads.out.bam
  } else {
    mark_duplicate_reads(merge_genome_bam.out.bam)
    bam_to_use = mark_duplicate_reads.out.bam
  }
  bam_to_fastq(bam_to_use)
  fastq = bam_to_fastq.out

  // Generate some useful files
  density(bam_to_use, star_ref_files)
  cufflinks(bam_to_use, annotation, sequins_iso_mix)
  stringtie(bam_to_use, annotation)
  feature_counts(bam_to_use, annotation)

  kallisto(fastq, kallisto_index, sequins_iso_mix)
  kallisto_advanced(fastq, kallisto_index, sequins_iso_mix)
  anaquin(bam_to_use, sequins_ref, kallisto_index, neat_mix_A, sequins_iso_mix)

  // QC Metrics
  insert_sizes(bam_to_use)
  rna_metrics(bam_to_use, flat_ref)
  adapter_count(bam_to_use)
  ribosomal_count(fastq)

}


process merge_transcriptome_bam {

  publishDir params.outdir, mode: "link"
  module "samtools/1.12"
  input:
    // Assume sorted by coord
    file("in*.bam")

  output:
    path("merged.transcriptome.bam")

  script:
    """
    numbam=\$(ls in*.bam | wc -l)
    if [ \$numbam -eq 1 ] ; then
      cp in*.bam "merged.transcriptome.bam"
    else
      samtools merge -n -f "merged.transcriptome.bam" in*.bam
    fi
    """
}

process merge_genome_bam {
  publishDir params.outdir, mode: "link"
  module "samtools/1.12"
  input:
    // Assume sorted by coord
    file("in*.bam")

  output:
    path("merged.genome.bam"), emit: bam
    path("*bai"), emit: bai

  script:
    """
    numbam=\$(ls in*.bam | wc -l)
    if [ \$numbam -eq 1 ] ; then
      cp in*.bam "merged.genome.bam"
      samtools index merged.genome.bam
    else
      samtools merge --write-index -f "merged.genome.bam##idx##merged.genome.bam.bai" in*.bam
    fi
    """
}


process remove_duplicate_reads {
  publishDir params.outdir, mode: "link"
  module "jdk/2.8.1", "picard/2.8.1", "samtools/1.12"
  label 'high_mem'
  input:
    path genomebam

  output:
    path "nodups.bam",                emit: bam
    path "picard.MarkDuplicates.txt", emit: metrics

  script:
    """
    # Add Mate Cigar information; required by UMI-aware MarkDuplicates
    # TODO: Actually required?
    picard RevertOriginalBaseQualitiesAndAddMateCigar \
      "INPUT=$genomebam" OUTPUT=cigar.bam \
      VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

    # Remove non-primary reads (also required)
    samtools view -F256 cigar.bam -o cigar_no_supp.bam

    picard UmiAwareMarkDuplicatesWithMateCigar \
      INPUT=cigar_no_supp.bam \
      OUTPUT=nodups.bam \
      METRICS_FILE=picard.MarkDuplicates.txt \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' \
      UMI_TAG_NAME=RX \
      REMOVE_DUPLICATES=true

    # Remove intermediate files
    rm cigar.bam cigar_no_supp.bam
  """

}

process mark_duplicate_reads {
  publishDir params.outdir, mode: "link"

  module "jdk/2.8.1", "picard/2.8.1", "samtools/1.12"
  input:
    path genomebam

  output:
    path "dupsmarked.bam",            emit: bam
    path "picard.MarkDuplicates.txt", emit: metrics

  script:
    """
    picard MarkDuplicates \
      "INPUT=${genomebam}" \
      OUTPUT=dupsmarked.bam \
      METRICS_FILE=picard.MarkDuplicates.txt \
      ASSUME_SORTED=true \
      VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    """

}


process bam_to_fastq {
  publishDir params.outdir, mode: "link"
  module "samtools/1.12"
  
  input:
    path input_bam
  output:
    tuple path("r1.fq.gz"), path("r2.fq.gz")

  script:
    """
    samtools sort -n -o namesorted.bam "$input_bam"
    samtools view -uf64 namesorted.bam \
      | samtools fastq - \
      | gzip -c \
      > r1.fq.gz
    samtools view -uf128 namesorted.bam \
      | samtools fastq - \
      | gzip -c \
      > r2.fq.gz
    """
}

process density {

  module "STAR/2.4.2a", "bedops/2.4.19", "htslib/1.6.0"
  publishDir params.outdir, mode: "link"
  input:
    path(input_bam)
    path("ref/*")

  output:
    path("Signal.*")

  script:
    """
    # Write starch and bigwig to .tmp files
    function convertBedGraph(){
      in="\$1"
      base="\$2"
      chrom="\$in.onlyChr.bg"
      grep '^chr' "\$in" | sort -k1,1 -k2,2n > \$chrom
      bedGraphToBigWig "\$chrom" chrNL.txt "\$base.bw"
      starch "\$chrom" > "\$base.starch"
    }

    mkdir -p Signal

    STAR --runMode inputAlignmentsFromBAM \
      --inputBAMfile "$input_bam" \
      --outWigType bedGraph \
      --outWigStrand Unstranded \
      --outFileNamePrefix Signal/ \
      --outWigReferencesPrefix chr \
      --outTmpDir STAR

    mv Signal/Signal.UniqueMultiple.str1.out.bg Signal/Signal.UniqueMultiple.unstranded.out.bg
    mv Signal/Signal.Unique.str1.out.bg         Signal/Signal.Unique.unstranded.out.bg

    STAR --runMode inputAlignmentsFromBAM \
      --inputBAMfile "$input_bam" \
      --outWigType bedGraph \
      --outWigStrand Stranded \
      --outFileNamePrefix Signal/ \
      --outWigReferencesPrefix chr \
      --outTmpDir STAR

    grep '^chr' ref/chrNameLength.txt > chrNL.txt

    convertBedGraph Signal/Signal.Unique.str1.out.bg               Signal.Unique.str-
    convertBedGraph Signal/Signal.Unique.str2.out.bg               Signal.Unique.str+
    convertBedGraph Signal/Signal.Unique.unstranded.out.bg         Signal.Unique.both
    convertBedGraph Signal/Signal.UniqueMultiple.str1.out.bg       Signal.UniqueMultiple.str-
    convertBedGraph Signal/Signal.UniqueMultiple.str2.out.bg       Signal.UniqueMultiple.str+
    convertBedGraph Signal/Signal.UniqueMultiple.unstranded.out.bg Signal.UniqueMultiple.both

    unstarch Signal.Unique.both.starch | bgzip > Signal.Unique.both.starch.bgz
    tabix -p bed Signal.Unique.both.starch.bgz
    """

}

process cufflinks {

  publishDir params.outdir, mode: "link"
  module "cufflinks/2.2.1", "R/3.2.5", "anaquin/2.0.1"

  input:
    path input_bam
    path annotation
    path sequins_iso_mix


  output:
    tuple path("genes.fpkm_tracking"), path("isoforms.fpkm_tracking"), path("anaquin_cufflinks/*")
    path "skipped.gtf"
    path "transcripts.gtf"


  script:
    """

    CUFF=cufflinks
    CUFF_COMMON="--no-update-check --library-type fr-firststrand"

    \$CUFF \$CUFF_COMMON --GTF "$annotation" "$input_bam"

    # delete duplicate rows and sort into identical orders across all samples
    Rscript \$STAMPIPES/scripts/rna-star/aggregate/dedupe_sort_cuffout.Rscript genes.fpkm_tracking
    Rscript \$STAMPIPES/scripts/rna-star/aggregate/dedupe_sort_cuffout.Rscript isoforms.fpkm_tracking
    mv genes.fpkm_tracking.sort genes.fpkm_tracking
    mv isoforms.fpkm_tracking.sort isoforms.fpkm_tracking

    # quantification with anaquin Rna Expression
    anaquin RnaExpression -o anaquin_cufflinks -rmix "$sequins_iso_mix" -usequin transcripts.gtf -mix A \
    || (echo "NA" > anaquin_cufflinks/RnaExpression_genes.tsv \
     && echo "NA" > anaquin_cufflinks/RnaExpression_isoforms.tsv \
     && echo "NA" > anaquin_cufflinks/RnaExpression_summary.stats)
    """

}

process stringtie {

  publishDir params.outdir, mode: "link"
  module "stringtie/1.3.4d"
  input:
    path(input_bam)
    path(annotation)

  output:
    path("stringtie_rf/*")
    

  script:
    """
    stringtie "$input_bam" --rf -p 4 -G "$annotation" \
      -o stringtie_rf/transcripts.gtf \
      -A stringtie_rf/abundances.txt
    """

}

process feature_counts {

    publishDir params.outdir, mode: "link"
    module "subread/1.5.1"

    input:
      path(input_bam)
      path(annotation)

    output:
      path "feature_counts.txt"
      path "feature_counts.txt.summary"

    script:
      """
      featureCounts --primary -B -C -p -P --fracOverlap .5 -s 2 -D 10000 \
        -t 'exon' -g 'gene_id' \
        -a "$annotation" -o feature_counts.txt "$input_bam"
      """
}


process kallisto {

  publishDir params.outdir, mode: "link"
  module "kallisto/0.43.1", "anaquin/2.0.1"

  input:
    tuple path(r1_fq), path(r2_fq)
    path kallisto_index
    path sequins_iso_mix

  output:
    path "anaquin_kallisto/*"
    path "kallisto_output/*"
    path "kallisto.log"

  script:
    """
    kallisto quant -i "${kallisto_index}" -o kallisto_output "${r1_fq}" "${r2_fq}" 2> kallisto.log

    anaquin RnaExpression -o anaquin_kallisto -rmix "${sequins_iso_mix}" -usequin kallisto_output/abundance.tsv -mix A \
    || (echo "NA" > anaquin_kallisto/RnaExpression_genes.tsv \
     && echo "NA" > anaquin_kallisto/RnaExpression_isoforms.tsv \
     && echo "NA" > anaquin_kallisto/RnaExpression_summary.stats)
    """
}

process kallisto_advanced {
  
  publishDir params.outdir, mode: "link"
  module "kallisto/0.43.1", "anaquin/2.0.1"

  input:
    tuple path(r1_fq), path(r2_fq)
    path kallisto_index
    path sequins_iso_mix

  output:
    path "anaquin_kallisto_adv/*"
    path "kallisto_output_adv/*"
    path "kallisto_adv.log"

  script:
    """
    kallisto quant --bias -b 100 --rf-stranded -i "${kallisto_index}" -o kallisto_output_adv "${r1_fq}" "${r2_fq}" 2> kallisto_adv.log

    anaquin RnaExpression -o anaquin_kallisto_adv -rmix "${sequins_iso_mix}" -usequin kallisto_output_adv/abundance.tsv -mix A \
    || (echo "NA" > anaquin_kallisto_adv/RnaExpression_genes.tsv \
     && echo "NA" > anaquin_kallisto_adv/RnaExpression_isoforms.tsv \
     && echo "NA" > anaquin_kallisto_adv/RnaExpression_summary.stats)
    """
}


process anaquin {

  module "samtools/1.12", "anaquin/2.0.1", "kallisto/0.43.1", "R/3.2.5"
  publishDir params.outdir, mode: "link"

  input:
    path input_bam
    path sequins_ref
    path kallisto_index
    path neat_mix_A
    path sequins_iso_mix

  output:
    path "anaquin_star/*"
    path "anaquin_subsample/anaquin.log"
    path "anaquin_subsample/RnaSubsample_summary.stats"
    path "anaquin_subsample/anaquin_kallisto/*"
    path "anaquin_subsample/anaquin_star/*"
    path "anaquin_subsample/kallisto_output/*"

  script:
  dilution = 0.0001
  """
  anaquin RnaAlign -rgtf "${sequins_ref}" -usequin "${input_bam}" -o anaquin_star
  bash \$STAMPIPES/scripts/rna-star/aggregate/anaquin_rnaalign_stats.bash anaquin_star/RnaAlign_summary.stats anaquin_star/RnaAlign_summary.stats.info

  # create subsample
  anaquin RnaSubsample -method "${dilution}" -usequin "${input_bam}" -o anaquin_subsample | samtools view -bS - > temp_subsample.bam

  # calculate anaquin align stats on subset alignment
  anaquin RnaAlign -rgtf "${sequins_ref}" -usequin temp_subsample.bam -o anaquin_subsample/anaquin_star
  bash \$STAMPIPES/scripts/rna-star/aggregate/anaquin_rnaalign_stats.bash \
    anaquin_subsample/anaquin_star/RnaAlign_summary.stats \
    anaquin_subsample/anaquin_star/RnaAlign_summary.stats.info

  # turn subset alignment to fastq
  samtools sort -n -o temp_subsample.sorted.bam temp_subsample.bam
  samtools view -uf64 temp_subsample.sorted.bam | samtools fastq - > subsample.fq1
  samtools view -uf128 temp_subsample.sorted.bam | samtools fastq - > subsample.fq2

  # call kallisto on subsampled fastqs
  kallisto quant -i "${kallisto_index}" -o anaquin_subsample/kallisto_output subsample.fq1 subsample.fq2

  # calculate kallisto stats on subsampled alignment
  anaquin RnaExpression -o anaquin_subsample/anaquin_kallisto -rmix "${sequins_iso_mix}" -usequin anaquin_subsample/kallisto_output/abundance.tsv -mix A \
    || (echo "NA" > anaquin_subsample/anaquin_kallisto/RnaExpression_genes.tsv \
     && echo "NA" > anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.tsv \
     && echo "NA" > anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats)

  bash $STAMPIPES/scripts/rna-star/aggregate/anaquin_rnaexp_stats.bash \
    anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats \
    anaquin_subsample/anaquin_kallisto/RnaExpression_summary.stats.info
  bash $STAMPIPES/scripts/rna-star/aggregate/anaquin_neat_comparison.bash \
    anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.tsv \
    anaquin_subsample/anaquin_kallisto/RnaExpression_isoforms.neatmix.tsv \
    "${neat_mix_A}"
"""

}

process insert_sizes {

  module "jdk", "picard/2.9.0", "R/3.2.5"
  publishDir params.outdir, mode: "link"

  input:
    path input_bam

  output:
    path "picard.CollectInsertSizes.txt"

  script:
    """
    picard CollectInsertSizeMetrics INPUT="${input_bam}" OUTPUT=picard.CollectInsertSizes.txt HISTOGRAM_FILE=/dev/null
    """
}

process rna_metrics {

  module "jdk", "picard/2.9.0"
  publishDir params.outdir, mode: "link"

  input:
    path input_bam
    path flat_ref

  output:
    path "picard.RnaSeqMetrics.txt"
    path "rna_stats_summary.info"

  script:
    """
    picard CollectRnaSeqMetrics INPUT="${input_bam}" OUTPUT=picard.RnaSeqMetrics.txt REF_FLAT="${flat_ref}" STRAND_SPECIFICITY="SECOND_READ_TRANSCRIPTION_STRAND"
    cat picard.RnaSeqMetrics.txt | grep -A 1 "METRICS CLASS" | sed 1d | tr '\t' '\n' > rna_stats_summary.info
    cat picard.RnaSeqMetrics.txt | grep -A 2 "METRICS CLASS" | sed 1d | sed 1d | tr '\t' '\n' | paste rna_stats_summary.info - > tmp.txt && mv tmp.txt rna_stats_summary.info
    """
}

process adapter_count {

  module "bwa", "samtools", "jdk", "picard/2.9.0"
  publishDir params.outdir, mode: "link"
  input:
    path input_bam

  output:
    path "adapter_counts.info"

  script:
    """
    adapter_count=\$(bash "\$STAMPIPES/scripts/bam/count_adapters.sh" "${input_bam}")
    if [ -n \$adapter_count ] ; then
            echo -e "adapter\t\$adapter_count" > "adapter_counts.info"
    fi
    """
}
process ribosomal_count {

  module "bowtie/1.0.0"
  publishDir params.outdir, mode: "link"

  input:
    tuple path(r1_fq), path(r2_fq)

  output:
    path("ribosomal_counts.info")

  script:
    """
    zcat -f "${r1_fq}" > trims.R1.fastq
    zcat -f "${r2_fq}" > trims.R2.fastq
    # TODO: Why is this hard-coded?
    num_reads=\$(bowtie -n 3 -e 140 /net/seq/data/genomes/human/GRCh38/noalts/contamination/hg_rRNA -1 trims.R1.fastq -2 trims.R2.fastq | wc -l )
    if [ -n \$num_reads ]; then
        echo -e "ribosomal-RNA\t\$num_reads" > "ribosomal_counts.info"
    fi
    """
}
