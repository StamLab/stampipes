/// Using the open-source adapter trimmer
/// https://github.com/OpenGene/fastp
process fastp_adapter_trim {

  // TODO: Get a newer version!
  // 0.23.0 enables deterministic results, which is crucial
  module 'fastp/0.21.0'
  cpus 3
  scratch false

  input:
    tuple path(r1), path(r2), val(adapterP5), val(adapterP7)

  output:
    path 'out.r?.fastq.gz', emit: fastq
    path 'fastp.json',      emit: metrics_json
    path 'fastp.html',      emit: metrics_html

  script:
    // TODO: Double-check adapter ordering
    """
    fastp \
      --in1 "${r1}" \
      --in2 "${r2}" \
      --adapter_sequence    "${adapterP7}" \
      --adapter_sequence_r2 "${adapterP5}" \
      --out1 "out.r1.fastq.gz" \
      --out2 "out.r2.fastq.gz" \
      --disable_quality_filtering \
      --disable_length_filtering \
      --thread 3
    """
}

process agent_adapter_trim {
  cpus 1
  scratch false
  module "jdk/11.0.16"

  input:
    tuple path(r1), path(r2)

  output:
    path 'output_R?.fastq.gz', emit: fastq

  script:
    """
    agent.sh trim -v2 \
      -fq1 "${r1}" \
      -fq2 "${r2}" \
      -out \$PWD/output
    """
}

/// Our custom in-house adapter-trimming script
process adapter_trim {
  cpus 3
  scratch false
  input:
    tuple path(r1), path(r2), val(adapterP5), val(adapterP7)

  output:
    path 'out.r?.fastq.gz', emit: fastq
    //path 'out.r2.fastq.gz', emit: trimmed_r2
    path 'adapter_trimming.txt', emit: metrics

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
