nextflow.enable.dsl=2

params.star_threads = 8
params.publish = true

process star {
  cpus params.star_threads
  module 'STAR/2.4.2a', 'samtools/1.7'

  publishDir params.outdir, enabled: params.publish

  label 'high_mem'
  //memory "256 GB"

  input:
    tuple path(r1_fq), path(r2_fq)
    path("ref/*")

  output:
    tuple path('Aligned.sortedByCoord.out.bam'), path('Aligned.toTranscriptome.out.bam')

  script:
    mode = "str_PE"
    threads = params.star_threads
    // Four-gigabyte fastq files will need more RAM to sort
    is_giant = r1_fq.size() > 3_000_000_000
    sort_ram = is_giant ? 200_000_000_000 : 30_000_000_000 
    sj_limit = is_giant ? "--limitOutSJcollapsed 5000000" : ""
    """
    # TODO: Update this??
    echo -e '@CO\tANNID:gencode.basic.tRNA.annotation.gtf.gz' > commentslong.txt

    STAR \
      --genomeDir "ref/"  \
      --readFilesIn "$r1_fq" "$r2_fq"   \
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
      --limitBAMsortRAM ${sort_ram} \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode TranscriptomeSAM \
      --outSAMheaderCommentFile commentslong.txt ${sj_limit}\
      --outSAMheaderHD '@HD' 'VN:1.4' 'SO:coordinate'
      # TODO: add stranded options
    """
}
