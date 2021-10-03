#!/usr/bin/env nextflow

params.r1 = ""
params.r2 = ""

params.adapter_p5 = ""
params.adapter_p7 = ""

params.seed = 12345
params.outdir = "output"

params.refdir = "/net/seq/data/genomes/human/GRCh38/noalts-sequins/"

starIndexDir = "${params.refdir}/STARgenome-gencode-v25"
rsemIndexDir = "${params.refdir}/RSEMgenome-gencode-v25"

params.star_threads = 8
params.rsem_threads = 8

process trim {

  input:
  file r1 from file(params.r1)
  file r2 from file(params.r2)
  val adapterP5 from params.adapter_p5
  val adapterP7 from params.adapter_p7

  output:
  set file('out.r1.fastq.gz'), file('out.r2.fastq.gz') into trimmed
  file('adapter_trimming.txt')

  script:
  """
  echo -e "P7\t$adapterP7\nP5\t$adapterP5" > adapters.txt

  trim-adapters-illumina \
    -f adapters.txt \
    -1 P5 -2 P7 \
    --threads=3 \
    "$r1" \
    "$r2" \
    out.r1.fastq.gz \
    out.r2.fastq.gz \
    &> adapter_trimming.txt
  """
}

process star {
  cpus params.star_threads
  module 'STAR', 'samtools/1.7'

  publishDir params.outdir, mode: "link"

  label 'high_mem'

  input:
  set file(r1), file(r2) from trimmed
  file 'ref/*' from Channel.fromPath("$starIndexDir/*").collect()
  val mode from "str_PE"
  val threads from params.star_threads

  output:
  file 'Aligned.sortedByCoord.out.bam' into coordinateBam
  file 'Aligned.toTranscriptome.out.bam' into transcriptBam

  script:

  """
  # TODO: Update this??
  echo -e '@CO\tANNID:gencode.basic.tRNA.annotation.gtf.gz' > commentslong.txt

  STAR \
    --genomeDir "ref/"  \
    --readFilesIn "$r1" "$r2"   \
    --outSAMunmapped Within --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD    \
    --outFilterMultimapNmax 20   \
    --outFilterMismatchNmax 999   \
    --outFilterMismatchNoverReadLmax 0.04   \
    --alignIntronMin 20   \
    --alignIntronMax 1000000   \
    --alignMatesGapMax 1000000   \
    --alignSJoverhangMin 8   \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --readFilesCommand zcat \
    --runThreadN ${threads} \
    --limitBAMsortRAM 30000000000 \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
    --outSAMheaderCommentFile commentslong.txt \
    --outSAMheaderHD '@HD' 'VN:1.4' 'SO:coordinate'
    # TODO: add stranded options

  """
}
