/// UMT-related processes

/// A generic script to move the UMT from a read name to the RX tag in a BAM file
process move_umt {
  input:
    tuple val(id), path(input_bam)
  output:
    tuple val(id), path("output.bam")
  shell:
    """
    python "$STAMPIPES/scripts/bam/move_umt_to_tag.py" \
      "$input_bam" \
      "output.bam"
    """
}

/// UMT-trimming for Takara Pico v3 kits
process takara_trim_umt {

  input:
    tuple path("in.r1.fq.gz"), path("in.r2.fq.gz")
    val readlength

  output:
    path "out.r*.fq.gz", emit: fastq
    path "takara_umt.log", emit: metrics
    
  script:
    """
    python "$STAMPIPES/scripts/fastq/takara_umt.py" \
      --readlength "$readlength" \
      <(zcat "in.r1.fq.gz") \
      <(zcat "in.r2.fq.gz") \
      >(gzip -c > "out.r1.fq.gz") \
      >(gzip -c > "out.r2.fq.gz") \
    > takara_umt.log

    """

}

