#!/usr/bin/env nextflow
/*
 * This is a proof-of-concept workflow for running our DNase alignment pipeline
 */

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

nuclear_chroms = "$params.genome" + ".nuclear.txt"
dataDir = "$baseDir/../../data"

def helpMessage() {
  log.info"""
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
    """.stripIndent();
}

if (params.help || !params.r1 || !params.r2 || !params.genome){
  helpMessage();
  exit 0;
}

// Some renaming for easier usage later
genome = params.genome
genome_name = file(params.genome).baseName
threads = params.threads

/*
 * Step 0: Split Fastq into chunks
 */
fastq_line_chunks = 4 * params.chunk_size
process split_r1_fastq {

  input:
  file(r1) from file(params.r1)
  val fastq_line_chunks

  output:
  file('split_r1*gz') into split_r1 mode flatten

  script:
  """
  zcat $r1 \
  | split -l $fastq_line_chunks \
    --filter='gzip -1 > \$FILE.gz' - 'split_r1'
  """
}
process split_r2_fastq {
  input:
  file(r2) from file(params.r2)
  val fastq_line_chunks

  output:
  file('split_r2*gz') into split_r2 mode flatten

  script:
  """
  zcat $r2 \
  | split -l $fastq_line_chunks \
    --filter='gzip -1 > \$FILE.gz' - 'split_r2'
  """
}

/*
 * Step 1: For each fastqc chunk, trim adapters
 */
process trim_adapters {

  cpus params.threads

  input:
  file split_r1
  file split_r2
  file adapters from file(params.adapter_file)

  output:
  set file('trim.R1.fastq.gz'), file('trim.R2.fastq.gz') into trimmed
  file('trim.counts.txt') into trim_counts

  script:
  """
  trim-adapters-illumina \
    -f "$adapters" \
    -1 P5 -2 P7 \
    --threads=${params.threads} \
    "$split_r1" \
    "$split_r2"  \
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

  input:
  set file(r1), file(r2) from trimmed

  output:
  set file('r1.trim.fastq.gz'), file('r2.trim.fastq.gz') into trimmed_fastq

  script:
  if (params.trim_to != 0)
    // TODO: Add padding to length with N's
    """
    zcat $r1 | awk 'NR%2==0 {print substr(\$0, 1, $params.trim_to)} NR%2!=0' | gzip -c -1 > r1.trim.fastq.gz
    zcat $r2 | awk 'NR%2==0 {print substr(\$0, 1, $params.trim_to)} NR%2!=0' | gzip -c -1 > r2.trim.fastq.gz
    """
  else
    """
    ln -s $r1 r1.trim.fastq.gz
    ln -s $r2 r2.trim.fastq.gz
    """

}

process add_umi_info {

  input:
  set file(r1), file(r2) from trimmed_fastq

  output:
  set file('r1.fastq.umi.gz'), file('r2.fastq.umi.gz') into with_umi

  script:
  if (params.UMI == 'thruplex')
    """
    python3 \$STAMPIPES/scripts/umi/extract_umt.py \
      <(zcat $r1) \
      <(zcat $r2) \
      >(gzip -c -1 > r1.fastq.umi.gz) \
      >(gzip -c -1 > r2.fastq.umi.gz)
    """
  else if (params.UMI == 'single-strand')
    """
    python3 \$STAMPIPES/scripts/umi/fastq_umi_add.py $r1 r1.fastq.umi.gz
    python3 \$STAMPIPES/scripts/umi/fastq_umi_add.py $r2 r2.fastq.umi.gz
    """

  else if (params.UMI == false || params.UMI == "")
    """
    ln -s $r1 r1.fastq.umi.gz
    ln -s $r2 r2.fastq.umi.gz
    """
  else
    error "--UMI must be `thruplex`, `single-strand` (for single-strand preparation), or false, got: '" + params.UMI + "'"
}

/*
 * Metrics: Fastq counts
 */
process fastq_counts {

  input:
  file(r1) from file(params.r1)
  file(r2) from file(params.r2)

  output:
  file 'fastq.counts' into fastq_counts

  script:
  """
  zcat $r1 \
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
  set file(trimmed_r1), file(trimmed_r2) from with_umi

  file genome from file(params.genome)
  file '*' from file("${params.genome}.amb")
  file '*' from file("${params.genome}.ann")
  file '*' from file("${params.genome}.bwt")
  file '*' from file("${params.genome}.fai")
  file '*' from file("${params.genome}.pac")
  file '*' from file("${params.genome}.sa")


  output:
  file 'out.bam' into unfiltered_bam

  script:
  """
  ls
  bwa aln \
    -Y -l 32 -n 0.04 \
    -t "${params.threads}" \
    "$genome" \
    "$trimmed_r1" \
    > out1.sai

  bwa aln \
    -Y -l 32 -n 0.04 \
    -t "$params.threads" \
    "$genome" \
    "$trimmed_r2" \
    > out2.sai

  bwa sampe \
    -n 10 -a 750 \
    "$genome" \
    out1.sai out2.sai \
    "$trimmed_r1" "$trimmed_r2" \
  | samtools view -b -t "$genome".fai - \
  > out.bam
  """

}

/*
 * Step 2b: filter bam files to have only good reads
 */
process filter_bam {

  input:
  file unfiltered_bam
  file nuclear_chroms from file(nuclear_chroms)

  output:
  file 'filtered.bam' into filtered_bam

  script:
  """
  python3 \$STAMPIPES/scripts/bwa/filter_reads.py \
  "$unfiltered_bam" \
  filtered.bam \
  "$nuclear_chroms"
  """
}

/*
 * Step 2c: sort bam files
 */
process sort_bam {

  cpus params.threads

  input:
  file filtered_bam

  output:
  file 'sorted.bam' into sorted_bam

  script:
  """
  samtools sort \
    -l 0 -m 1G -@ "${params.threads}" "$filtered_bam" \
    > sorted.bam
  """
}

/*
 * Step 3: Merge alignments into one big ol' file
 */
process merge_bam {
  input:
  file 'sorted_bam_*' from sorted_bam.collect()

  output:
  file 'merged.bam' into merged_bam

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

  publishDir params.outdir, mode: "link"

  input:
  file(merged_bam) from merged_bam

  output:
  file 'marked.bam' into marked_bam
  file 'marked.bam' into marked_bam_for_counts
  file 'MarkDuplicates.picard'


  script:
  if (params.UMI)
    """
    picard RevertOriginalBaseQualitiesAndAddMateCigar \
      INPUT=$merged_bam OUTPUT=cigar.bam \
      VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

    picard UmiAwareMarkDuplicatesWithMateCigar INPUT=cigar.bam OUTPUT=marked.bam \
      METRICS_FILE=MarkDuplicates.picard UMI_TAG_NAME=XD ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    """
  else
    """
    picard RevertOriginalBaseQualitiesAndAddMateCigar \
      INPUT=$merged_bam OUTPUT=cigar.bam \
      VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

    picard MarkDuplicatesWithMateCigar INPUT=cigar.bam OUTPUT=marked.bam \
      METRICS_FILE=MarkDuplicates.picard ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' \
      MINIMUM_DISTANCE=300
    """
}

/*
 * Step 5: Filter bam file
 */
filter_flag = 512
if (params.UMI)
  filter_flag = 1536

process filter_bam_to_unique {

  publishDir params.outdir, mode: "link"

  input:
  file marked_bam

  output:
  set file('filtered.bam'), file('filtered.bam.bai') into uniquely_mapping_bam

  script:
  """
  samtools view $marked_bam -b -F $filter_flag > filtered.bam
  samtools index filtered.bam
  """

}

// Marked bam file is used by several downstream results
uniquely_mapping_bam.into { bam_for_insert; bam_for_spot; bam_for_density }

/*
 * Metrics: bam counts
 */
process bam_counts {


  input:
  file(sorted_bam) from marked_bam_for_counts

  output:
  file('bam.counts.txt') into bam_counts

  script:
  """
  python3 \$STAMPIPES/scripts/bwa/bamcounts.py \
    "$sorted_bam" \
    bam.counts.txt
  """
}

/*
 * Metrics: Insert size
 */
process insert_size {

  publishDir params.outdir, mode: "link"

  input:
  set file(bam), file(bai) from bam_for_insert

  output:
  file 'CollectInsertSizeMetrics.picard'
  file 'CollectInsertSizeMetrics.picard.pdf'

  script:
  """
  samtools idxstats "$bam" \
  | cut -f 1 \
  | grep -v chrM \
  | grep -v chrC \
  | xargs samtools view -b "$bam" \
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

  publishDir params.outdir, mode: "link"

  input:
  set file(bam), file(bai) from bam_for_spot
  file('*') from file("${dataDir}/annotations/${genome_name}.K${params.readlength}.mappable_only.bed")
  file('*') from file("${dataDir}/annotations/${genome_name}.chromInfo.bed")

  output:
  file 'subsample.r1.spot.out'
  file 'spotdups.txt'

  script:
  """
  # random sample
  samtools view -h -F 12 -f 3 "$bam" \
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
win = 75
bini = 20
process density_files {

  publishDir params.outdir, mode: "link"

  input:
  set file(bam), file(bai) from bam_for_density
  file fai from file("${params.genome}.fai")
  file density_buckets from file(
    "$baseDir/../../data/densities/chrom-buckets.${genome_name}.${win}_${bini}.bed.starch"
  )

  output:
  file 'density.bed.starch'
  file 'density.bw'
  file 'density.bed.bgz'


  script:
  """
  bam2bed -d \
    < $bam \
    | cut -f1-6 \
    | awk '{ if( \$6=="+" ){ s=\$2; e=\$2+1 } else { s=\$3-1; e=\$3 } print \$1 "\t" s "\t" e "\tid\t" 1 }' \
    | sort-bed - \
    > sample.bed

  unstarch "${density_buckets}" \
    | bedmap --faster --echo --count --delim "\t" - sample.bed \
    | awk -v binI=$bini -v win=$win \
        'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print \$1 "\t" \$2 + shiftFactor "\t" \$3-shiftFactor "\tid\t" i \$4}' \
    | starch - \
    > density.bed.starch

  unstarch density.bed.starch | awk -v binI=$bini -f \$STAMPIPES/awk/bedToWig.awk > density.wig

  wigToBigWig -clip density.wig "${fai}" density.bw


  unstarch density.bed.starch | bgzip > density.bed.bgz
  tabix -p bed density.bed.bgz
  """
}

/*
 * Metrics: total counts
 */
process total_counts {

  publishDir params.outdir, mode: "link"

  input:
  file 'fastqcounts*' from fastq_counts.collect()
  file 'trimcounts*' from trim_counts.collect()
  file 'bamcounts*' from bam_counts.collect()

  output:
  file 'tagcounts.txt'

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
