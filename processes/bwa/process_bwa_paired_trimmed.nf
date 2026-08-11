#!/usr/bin/env nextflow
/*
 * This is our main DNase alignment pipeline
 */

nextflow.enable.dsl = 2


params.help = false
params.threads = 1
params.chunk_size = 16000000

params.UMI = false
params.trim_to = 0
params.genome = ""
params.r1 = ""
params.r2 = ""
params.outdir = "output"
params.readlength = 36

params.cramthreads = 10

def helpMessage() {
  log.info(
    """
    Usage: nextflow run process_bwa_paired_trimmed.nf \\
             --r1 r1.fastq.gz \\
             --r2 r2.fastq.gz \\
             --adapter_file adapters.txt \\
             --genome /path/to/genome \\
    Options:
    --threads [count]     The number of threads that will be used for child applications  (1)
    --chunk_size [count]  How many reads to process per chunk                             (16000000)
    --UMI                 The reads contain UMI markers ('single-strand', 'thruplex')         (false)
    --trim_to [length]    Trim fastq reads to [length] bp (0 for no trimming)             (0)
    --outdir [dir]        Where to write results to                                       (output)
    """.stripIndent()
  )
}


process split_r1_fastq {
  input:
  path r1
  val fastq_line_chunks

  output:
  path 'split_r1*gz'

  script:
  """
  zcat ${r1} \
  | split -l ${fastq_line_chunks} \
    --filter='gzip -1 > \$FILE.gz' - 'split_r1'
  """
}
process split_r2_fastq {
  input:
  path r2
  val fastq_line_chunks

  output:
  path 'split_r2*gz'

  script:
  """
  zcat ${r2} \
  | split -l ${fastq_line_chunks} \
    --filter='gzip -1 > \$FILE.gz' - 'split_r2'
  """
}

/*
 * Step 1: For each fastqc chunk, trim adapters
 */
process trim_adapters {

  cpus params.threads
  scratch false

  input:
  path split_r1
  path split_r2
  path adapters

  output:
  tuple path('trim.R1.fastq.gz'), path('trim.R2.fastq.gz')
  path 'trim.counts.txt'

  script:
  """
  trim-adapters-illumina \
    -f "${adapters}" \
    -1 P5 -2 P7 \
    --threads=${params.threads} \
    "${split_r1}" \
    "${split_r2}"  \
    "trim.R1.fastq.gz" \
    "trim.R2.fastq.gz" \
  &> trimstats.txt

  awk '{print "adapter-trimmed\t" \$NF * 2}' \
  < trimstats.txt \
  > trim.counts.txt
  """
}

/*
 * Step 1.1: Trim to the appropriate length
 */
process trim_to_length {

  scratch false

  input:
  tuple path(r1), path(r2)

  output:
  tuple path('r1.trim.fastq.gz'), path('r2.trim.fastq.gz')

  script:
  if (params.trim_to != 0) {
    // TODO: Add padding to length with N's
    """
    zcat ${r1} | awk 'NR%2==0 {print substr(\$0, 1, ${params.trim_to})} NR%2!=0' | gzip -c -1 > r1.trim.fastq.gz
    zcat ${r2} | awk 'NR%2==0 {print substr(\$0, 1, ${params.trim_to})} NR%2!=0' | gzip -c -1 > r2.trim.fastq.gz
    """
  }
  else {
    """
    ln -s ${r1} r1.trim.fastq.gz
    ln -s ${r2} r2.trim.fastq.gz
    """
  }
}

process add_umi_info {

  scratch false

  input:
  tuple path(r1), path(r2)

  output:
  tuple path('r1.fastq.umi.gz'), path('r2.fastq.umi.gz')

  script:
  if (params.UMI == 'thruplex') {
    """
    python3 \$STAMPIPES/scripts/umi/extract_umt.py \
      <(zcat ${r1}) \
      <(zcat ${r2}) \
      >(gzip -c -1 > r1.fastq.umi.gz) \
      >(gzip -c -1 > r2.fastq.umi.gz)
    """
  }
  else if (params.UMI == 'single-strand') {
    """
    python3 \$STAMPIPES/scripts/umi/fastq_umi_add.py ${r1} r1.fastq.umi.gz
    python3 \$STAMPIPES/scripts/umi/fastq_umi_add.py ${r2} r2.fastq.umi.gz
    """
  }
  else if (params.UMI == false || params.UMI == "") {
    """
    ln -s ${r1} r1.fastq.umi.gz
    ln -s ${r2} r2.fastq.umi.gz
    """
  }
  else {
    error("--UMI must be `thruplex`, `single-strand` (for single-strand preparation), or false, got: '" + params.UMI + "'")
  }
}

/*
 * Metrics: Fastq counts
 */
process fastq_counts {

  scratch false

  input:
  path r1
  path r2

  output:
  path 'fastq.counts'

  script:
  """
  zcat ${r1} \
  | awk -v paired=1 -f \$STAMPIPES/awk/illumina_fastq_count.awk \
  > fastq.counts
  """
}

/*
 * Step 2a: Create alignment files
 */
process align {

  cpus params.threads

  input:
  tuple path(trimmed_r1), path(trimmed_r2)
  path genome
  path '*'
  path '*'
  path '*'
  path '*'
  path '*'
  path '*'

  output:
  path 'out.bam'

  script:
  """
  ls
  bwa aln \
    -Y -l 32 -n 0.04 \
    -t "${params.threads}" \
    "${genome}" \
    "${trimmed_r1}" \
    > out1.sai

  bwa aln \
    -Y -l 32 -n 0.04 \
    -t "${params.threads}" \
    "${genome}" \
    "${trimmed_r2}" \
    > out2.sai

  bwa sampe \
    -n 10 -a 750 \
    "${genome}" \
    out1.sai out2.sai \
    "${trimmed_r1}" "${trimmed_r2}" \
  | samtools view -b -t "${genome}".fai - \
  > out.bam
  """
}

/*
 * Step 2b: filter bam files to have only good reads
 */
process filter_bam {

  scratch false

  input:
  path unfiltered_bam
  path nuclear_chroms

  output:
  path 'filtered.bam'

  script:
  """
  python3 \$STAMPIPES/scripts/bwa/filter_reads.py \
  "${unfiltered_bam}" \
  filtered.bam \
  "${nuclear_chroms}"
  """
}

/*
 * Step 2c: sort bam files
 */
process sort_bam {

  cpus params.threads
  scratch false

  input:
  path filtered_bam

  output:
  path 'sorted.bam'

  script:
  """
  samtools sort \
    -l 0 -m 2G -@ "${params.threads}" "${filtered_bam}" \
    > sorted.bam
  """
}

/*
 * Step 3: Merge alignments into one big ol' file
 */
process merge_bam {
  scratch false

  input:
  path 'sorted_bam_*'

  output:
  path 'merged.bam'

  script:
  """
  samtools merge merged.bam sorted_bam*
  samtools index merged.bam
  """
}

/*
 * Step 4: Mark duplicates with Picard
 */
process mark_duplicates {

  label "high_mem"

  publishDir params.outdir

  input:
  path merged_bam

  output:
  path 'marked.bam'
  path 'MarkDuplicates.picard'

  script:
  if (params.UMI) {
    """
    picard RevertOriginalBaseQualitiesAndAddMateCigar \
      INPUT=${merged_bam} OUTPUT=cigar.bam \
      VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

    picard UmiAwareMarkDuplicatesWithMateCigar INPUT=cigar.bam OUTPUT=marked.bam \
      METRICS_FILE=MarkDuplicates.picard UMI_TAG_NAME=XD ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    """
  } else {
    """
    picard RevertOriginalBaseQualitiesAndAddMateCigar \
      INPUT=${merged_bam} OUTPUT=cigar.bam \
      VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

    picard MarkDuplicatesWithMateCigar INPUT=cigar.bam OUTPUT=marked.bam \
      METRICS_FILE=MarkDuplicates.picard ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' \
      MINIMUM_DISTANCE=300
    """
  }
}

/*
 * Step 5: Filter bam file
 */


process filter_bam_to_unique {

  scratch false

  input:
  path marked_bam

  output:
  tuple path('filtered.bam'), path('filtered.bam.bai')

  script:
  filter_flag = params.UMI ? 1536 : 512
  """
  samtools view ${marked_bam} -b -F ${filter_flag} > filtered.bam
  samtools index filtered.bam
  """
}

/*
 * Metrics: bam counts
 */
process bam_counts {

  scratch false

  input:
  path sorted_bam

  output:
  path 'bam.counts.txt'

  script:
  """
  python3 \$STAMPIPES/scripts/bwa/bamcounts.py \
    "${sorted_bam}" \
    bam.counts.txt
  """
}

/*
 * Metrics: Insert size
 */
process insert_size {

  publishDir params.outdir

  input:
  tuple path(bam), path(bai)

  output:
  path 'CollectInsertSizeMetrics.picard'
  path 'CollectInsertSizeMetrics.picard.pdf'

  script:
  """
  samtools idxstats "${bam}" \
  | cut -f 1 \
  | grep -v chrM \
  | grep -v chrC \
  | xargs samtools view -b "${bam}" \
  > nuclear.bam

  picard CollectInsertSizeMetrics \
    INPUT=nuclear.bam \
    OUTPUT=CollectInsertSizeMetrics.picard \
    HISTOGRAM_FILE=CollectInsertSizeMetrics.picard.pdf \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true
  """
}

/*
 * Metrics: SPOT score
 */

process spot_score {

  publishDir params.outdir

  input:
  tuple path(bam), path(bai)
  path '*'
  path '*'

  output:
  path 'subsample.r1.spot.out'
  path 'spotdups.txt'

  script:
  genome_name = file(params.genome).baseName
  """
  # random sample
  samtools view -h -F 12 -f 3 "${bam}" \
    | awk '{if( ! index(\$3, "chrM") && \$3 != "chrC" && \$3 != "random"){print}}' \
    | samtools view -1 - \
    -o paired.bam
  bash \$STAMPIPES/scripts/bam/random_sample.sh paired.bam subsample.bam 5000000
        samtools view -1 -f 0x0040 subsample.bam -o subsample.r1.bam

  # hotspot
  bash \$STAMPIPES/scripts/SPOT/runhotspot.bash \
    \$HOTSPOT_DIR \
    \$PWD \
    \$PWD/subsample.r1.bam \
    "${genome_name}" \
    "${params.readlength}" \
    DNaseI

  # Remove existing duplication marks
  picard RevertSam \
    INPUT=subsample.bam \
    OUTPUT=clear.bam \
    VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATE_INFORMATION=true SORT_ORDER=coordinate \
    RESTORE_ORIGINAL_QUALITIES=false REMOVE_ALIGNMENT_INFORMATION=false

	picard MarkDuplicatesWithMateCigar \
    INPUT=clear.bam \
    METRICS_FILE=spotdups.txt \
    OUTPUT=/dev/null \
    ASSUME_SORTED=true \
    MINIMUM_DISTANCE=300 \
    VALIDATION_STRINGENCY=SILENT \
		READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
  """
}

/*
 * Density tracks (for browser)
 */
process density_files {

  label "high_mem"
  publishDir params.outdir

  input:
  tuple path(bam), path(bai)
  path fai
  path density_buckets

  output:
  path 'density.bed.starch'
  path 'density.bw'
  path 'density.bed.bgz'

  script:
  win = 75
  bini = 20
  """
  ls -l
  bam2bed -d \
    < ${bam} \
    | cut -f1-6 \
    | awk '{ if( \$6=="+" ){ s=\$2; e=\$2+1 } else { s=\$3-1; e=\$3 } print \$1 "\t" s "\t" e "\tid\t" 1 }' \
    | sort-bed - \
    > sample.bed

  unstarch "${density_buckets}" \
    | bedmap --faster --echo --count --delim "\t" - sample.bed \
    | awk -v binI=${bini} -v win=${win} \
        'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print \$1 "\t" \$2 + shiftFactor "\t" \$3-shiftFactor "\tid\t" i \$4}' \
    | starch - \
    > density.bed.starch

  unstarch density.bed.starch | awk -v binI=${bini} -f \$STAMPIPES/awk/bedToWig.awk > density.wig
  echo "########"
  ls -l
  wigToBigWig -clip density.wig "${fai}" density.bw


  unstarch density.bed.starch | bgzip > density.bed.bgz
  tabix -p bed density.bed.bgz
  """
}

/*
 * Metrics: total counts
 */
process total_counts {

  publishDir params.outdir

  input:
  path 'fastqcounts*'
  path 'trimcounts*'
  path 'bamcounts*'

  output:
  path 'tagcounts.txt'

  script:
  """
  cat *counts* \
  | awk '
  { x[\$1] += \$2 }
  END {for (i in x) print i "\t" x[i]}
  ' \
  | sort -k 1,1 \
  > tagcounts.txt
  """
}

process cram {
  publishDir params.outdir
  cpus params.cramthreads / 2
  scratch false

  input:
  path bam
  path ref
  path fai

  output:
  path cramfile
  path "${cramfile}.crai"

  script:
  cramfile = bam.name.replace("bam", "cram")
  """
  samtools view "${bam}" \
    -C -O cram,version=3.0,level=7,lossy_names=0 \
    -T "${ref}" \
    --threads "${params.cramthreads}" \
    --write-index \
    -o "${cramfile}"
  """
}

workflow {
  // Check for required parameters
  if (params.help || !params.r1 || !params.r2 || !params.genome) {
    helpMessage()
    exit(0)
  }

  // Setup variables
  nuclear_chroms = "${params.genome}.nuclear.txt"
  dataDir = "${baseDir}/../../data"
  genome_name = file(params.genome).baseName
  fastq_line_chunks = 4 * params.chunk_size

  // Input channels
  r1_ch = Channel.fromPath(params.r1)
  r2_ch = Channel.fromPath(params.r2)
  adapter_file_ch = Channel.fromPath(params.adapter_file)

  // Genome files channels
  genome_ch = Channel.fromPath(params.genome)
  genome_amb_ch = Channel.fromPath("${params.genome}.amb")
  genome_ann_ch = Channel.fromPath("${params.genome}.ann")
  genome_bwt_ch = Channel.fromPath("${params.genome}.bwt")
  genome_fai_ch = Channel.fromPath("${params.genome}.fai")
  genome_pac_ch = Channel.fromPath("${params.genome}.pac")
  genome_sa_ch = Channel.fromPath("${params.genome}.sa")
  genome_fa_ch = Channel.fromPath("${params.genome}.fa")
  nuclear_chroms_ch = Channel.fromPath(nuclear_chroms)

  // Annotation files
  mappable_bed_ch = Channel.fromPath("${dataDir}/annotations/${genome_name}.K${params.readlength}.mappable_only.bed")
  chrominfo_bed_ch = Channel.fromPath("${dataDir}/annotations/${genome_name}.chromInfo.bed")
  density_buckets_ch = Channel.fromPath("${baseDir}/../../data/densities/chrom-buckets.${genome_name}.75_20.bed.starch")

  // Step 0: Split fastq files
  split_r1_fastq(r1_ch, fastq_line_chunks)
  split_r2_fastq(r2_ch, fastq_line_chunks)

  // Step 1: Trim adapters 
  trimmed_reads = trim_adapters(
    split_r1_fastq.out.flatten(),
    split_r2_fastq.out.flatten(),
    adapter_file_ch,
  )

  // Step 1.1: Trim to length
  length_trimmed = trim_to_length(trimmed_reads[0])

  // Step 1.2: Add UMI info
  umi_reads = add_umi_info(length_trimmed)

  // Metrics: Fastq counts
  fastq_count_results = fastq_counts(r1_ch, r2_ch)

  // Step 2a: Align reads
  unfiltered_bams = align(
    umi_reads,
    genome_ch,
    genome_amb_ch,
    genome_ann_ch,
    genome_bwt_ch,
    genome_fai_ch,
    genome_pac_ch,
    genome_sa_ch,
  )

  // Step 2b: Filter bam files
  filtered_bams = filter_bam(unfiltered_bams, nuclear_chroms_ch)

  // Step 2c: Sort bam files
  sorted_bams = sort_bam(filtered_bams)

  // Step 3: Merge alignments
  merged_bam = merge_bam(sorted_bams.collect())

  // Step 4: Mark duplicates
  mark_dup_results = mark_duplicates(merged_bam)
  marked_bam = mark_dup_results[0]

  // Step 5: Filter to unique reads
  unique_bam = filter_bam_to_unique(marked_bam)

  // Metrics: BAM counts
  bam_count_results = bam_counts(marked_bam)

  // Metrics: Insert size
  insert_size(unique_bam)

  // Metrics: SPOT score
  spot_score(unique_bam, mappable_bed_ch, chrominfo_bed_ch)

  // Density files
  density_files(unique_bam, genome_fai_ch, density_buckets_ch)

  // Total counts
  total_counts(
    fastq_count_results.collect(),
    trimmed_reads[1].collect(),
    bam_count_results.collect(),
  )

  // CRAM conversion
  cram(marked_bam, genome_fa_ch, genome_fai_ch)
}
