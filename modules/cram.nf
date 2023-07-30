nextflow.enable.dsl=2
params.cram_compression_threads = 10
params.cram_version = "3.0"
params.cram_compression_level = 7
params.cram_lossy_names = 0
params.cram_write_index = true
params.cram_compress_other_args = ""  // format:  "key1=value,key2=value2", etc.

workflow test {
  def meta = [id: "test"]
  def bam = file("test_data/dnase/aggregation/in_1.bam" )
  def ref = file("test_data/ref/chr22.fa")
  encode_cram([meta, bam, ref])
  encode_cram_no_ref([meta, bam])

  encode_cram.out.cram.view { "cram: ${it}" }
  encode_cram.out.cram_with_index.view { "cram with index: ${it}" }
}

process encode_cram {

  tag "${meta?.id}"

  module "samtools/1.12"
  container "quay.io/biocontainers/samtools:1.12--h9aed4be_1"

  cpus Math.ceil(params.cram_compression_threads / 2)
  scratch false

  input:
    tuple val(meta), path(input_bam), path(reference)

  output:
    tuple val(meta), path(output_cram_name), emit: cram
    tuple val(meta), path(output_cram_name), path("${output_cram_name}.crai"), emit: cram_with_index optional true

  script:
    output_cram_name = "${input_bam.baseName}.cram"
    fmt_options = [
      "version=${params.cram_version}",
      "level=${params.cram_compression_level}",
      "lossy_names=${params.cram_lossy_names}",
      "${params.cram_compress_other_args}",
      ].join(",")
    writeindexflag = params.cram_write_index ? "--write-index" : ""
    """
    samtools view \
      --reference "${reference}" \
      -C -o "${output_cram_name}" \
      --output-fmt-option "${fmt_options}" \
      --threads "${params.cram_compression_threads}" \
      ${writeindexflag} \
      "${input_bam}"
    """
}

process sort_and_encode_cram_no_ref {

  tag "${meta.id}"

  module "samtools/1.12"
  container "quay.io/biocontainers/samtools:1.12--h9aed4be_1"

  cpus Math.ceil(params.cram_compression_threads / 2)

  input:
    tuple val(meta), path(input_bam)

  output:
    tuple val(meta), path(output_cram_name), emit: cram
    tuple val(meta), path(output_cram_name), path("${output_cram_name}.crai"), emit: cram_with_index optional true

  script:
    output_cram_name = "${input_bam.baseName}.cram"
    fmt_options = [
      "no_ref=1",
      "version=${params.cram_version}",
      "level=${params.cram_compression_level}",
      "lossy_names=${params.cram_lossy_names}",
      "${params.cram_compress_other_args}",
      ].join(",")
    writeindexflag = params.cram_write_index ? "--write-index" : ""
    """
    samtools sort \
      -o "${output_cram_name}" \
      --output-fmt-option "${fmt_options}" \
      --threads "${params.cram_compression_threads}" \
      ${writeindexflag} \
      "${input_bam}"
    """
}


process sort_and_encode_cram {

  tag "${meta.id}"

  module "samtools/1.12"
  container "quay.io/biocontainers/samtools:1.12--h9aed4be_1"

  cpus Math.ceil(params.cram_compression_threads / 2)

  input:
    tuple val(meta), path(input_bam), path(reference)

  output:
    tuple val(meta), path(output_cram_name), emit: cram
    tuple val(meta), path(output_cram_name), path("${output_cram_name}.crai"), emit: cram_with_index optional true

  script:
    output_cram_name = "${input_bam.baseName}.cram"
    fmt_options = [
      "version=${params.cram_version}",
      "level=${params.cram_compression_level}",
      "lossy_names=${params.cram_lossy_names}",
      "${params.cram_compress_other_args}",
      ].join(",")
    writeindexflag = params.cram_write_index ? "--write-index" : ""

    """
    samtools sort \
      --reference "${reference}" \
      -o "${output_cram_name}" \
      --output-fmt-option "${fmt_options}" \
      --threads "${params.cram_compression_threads}" \
      ${writeindexflag} \
      "${input_bam}"
    """
}

process encode_cram_no_ref {
  tag "${meta.id}"

  module "samtools/1.12"
  container "quay.io/biocontainers/samtools:1.12--h9aed4be_1"

  cpus Math.ceil(params.cram_compression_threads / 2)
  scratch false

  input:
    tuple val(meta), path(input_bam)

  output:
    tuple val(meta), path(output_cram_name), emit: cram
    tuple val(meta), path(output_cram_name), path("${output_cram_name}.crai"), emit: cram_with_index optional true

  script:
    output_cram_name = "${input_bam.baseName}.cram"
    fmt_options = [
      "no_ref=1",
      "version=${params.cram_version}",
      "level=${params.cram_compression_level}",
      "lossy_names=${params.cram_lossy_names}",
      "${params.cram_compress_other_args}",
      ].join(",")
    writeindexflag = params.cram_write_index ? "--write-index" : ""
    """
    samtools view \
      -C -o "${output_cram_name}" \
      --output-fmt-option "${fmt_options}" \
      --threads "${params.cram_compression_threads}" \
      ${writeindexflag} \
      "${input_bam}"
    """
}
