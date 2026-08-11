nextflow.enable.dsl = 2

params.help = false
params.threads = 1

params.UMI = false
params.genome = ""
params.outdir = "output"
params.domotifs = false
params.dofeatures = false

params.readlength = 36

params.peakcaller = "hotspot2"

params.bias = ""
params.chunksize = 5000

params.hotspot_id = "default"
params.hotspot_index = "."

params.cramthreads = 10

def helpMessage() {
  log.info(
    """
    Usage: nextflow run basic.nf \\
             --genome /path/to/genome \\
             --bams '1.bam,2.bam...' \\
             --UMI true/false        \\
             --outdir /path/to/output

    """.stripIndent()
  )
}


process merge_bams {
  label "modules"

  input:
  // Inputs may be cram file but that's fine
  path 'in*.bam'

  output:
  path 'merged.bam'

  script:
  """
  samtools merge 'merged.bam' in*.bam
  """
}

// TODO: single end
process dups {
  label "modules"
  label 'high_mem'
  publishDir params.outdir, saveAs: { name -> name ==~ /.*ba[mi]$/ ? null : name }

  input:
  path merged

  output:
  path 'marked.bam', emit: marked_bam
  path 'marked.bam', emit: bams_to_cram_marked
  path 'marked.bam.bai'
  path 'MarkDuplicates.picard'

  script:
  if (params.UMI) {
    cmd = "UmiAwareMarkDuplicatesWithMateCigar"
    extra = "UMI_TAG_NAME=XD"
  }
  else {
    cmd = "MarkDuplicatesWithMateCigar"
    extra = "MINIMUM_DISTANCE=300"
  }
  """
  env
  ls -l
  picard RevertOriginalBaseQualitiesAndAddMateCigar \
    "INPUT=${merged}" OUTPUT=cigar.bam \
    VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0
  echo "###"
  ls -l
  picard "${cmd}" \
      INPUT=cigar.bam OUTPUT=marked.bam \
      ${extra} \
      METRICS_FILE=MarkDuplicates.picard ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
  echo "###"
  ls -l

  samtools index marked.bam
  """
}

process filter_bam {
  label "modules"

  input:
  path bam

  output:
  path "filtered.bam", emit: filtered_bam

  script:
  flag = params.UMI ? 1536 : 512
  """
  samtools view -b -F "${flag}" marked.bam > filtered.bam
  """
}

process filter_nuclear {
  label "modules"

  input:
  path bam
  path nuclear_chroms

  output:
  path 'nuclear.bam', emit: nuclear_bam

  script:
  """
  samtools index "${bam}"
  cat "${nuclear_chroms}" \
  | xargs samtools view -b "${bam}" \
  > nuclear.bam
  """
}

process macs2 {
  label "macs2"
  publishDir "${params.outdir}/peaks_macs2"
  scratch false

  input:
  path bam

  output:
  path 'NA_*'

  when:
  params.peakcaller == "macs2"

  script:
  """
  macs2 callpeak \
    -t "${bam}" \
    -f BAMPE \
    -g hs \
    -B \
    -q 0.01
  """
}

process hotspot2 {
  label "modules"
  label 'high_mem'
  publishDir "${params.outdir}"
  container "fwip/hotspot2:latest"

  input:
  val hotspotid
  path nuclear
  path mappable
  path chrom_sizes
  path centers

  output:
  path 'peaks/nuclear*'
  path 'peaks/nuclear.hotspots.fdr0.05.starch', emit: hotspot_calls
  path 'peaks/nuclear.hotspots.fdr0.05.starch', emit: hotspot_calls_for_bias
  path 'peaks/nuclear.peaks.fdr0.001.starch', emit: onepercent_peaks

  when:
  params.peakcaller == "hotspot2"

  script:
  """
  export TMPDIR=\$(mktemp -d)
  hotspot2.sh -F 0.5 -p "varWidth_20_${hotspotid}" \
    -M "${mappable}" \
    -c "${chrom_sizes}" \
    -C "${centers}" \
    "${nuclear}" \
    'peaks'

  cd peaks

  # Rename peaks files to include FDR
  mv nuclear.peaks.narrowpeaks.starch nuclear.peaks.narrowpeaks.fdr0.05.starch
  mv nuclear.peaks.starch nuclear.peaks.fdr0.05.starch

  bash \$STAMPIPES/scripts/SPOT/info.sh \
    nuclear.hotspots.fdr0.05.starch hotspot2 nuclear.SPOT.txt \
    > nuclear.hotspot2.info

  # TODO: Move this to separate process
  hsmerge.sh -f 0.01 nuclear.allcalls.starch nuclear.hotspots.fdr0.01.starch
  hsmerge.sh -f 0.001 nuclear.allcalls.starch nuclear.hotspots.fdr0.001.starch

  density-peaks.bash \$TMPDIR "varWidth_20_${hotspotid}" nuclear.cutcounts.starch nuclear.hotspots.fdr0.01.starch ../"${chrom_sizes}" nuclear.density.starch nuclear.peaks.fdr0.01.starch \$(cat nuclear.cleavage.total)
  density-peaks.bash \$TMPDIR "varWidth_20_${hotspotid}" nuclear.cutcounts.starch nuclear.hotspots.fdr0.001.starch ../"${chrom_sizes}" nuclear.density.starch nuclear.peaks.fdr0.001.starch \$(cat nuclear.cleavage.total)

  rm -rf "\$TMPDIR"
  """
}

process spot_score {
  label "modules"
  publishDir params.outdir

  input:
  path bam
  path mappable
  path chromInfo
  val genome_name

  output:
  path 'r1.spot.out'
  path 'r1.hotspot.info'

  script:
  """
  # random sample
	samtools view -h -F 12 -f 3 "${bam}" \
		| awk '{if( ! index(\$3, "chrM") && \$3 != "chrC" && \$3 != "random"){print}}' \
		| samtools view -uS - \
		> nuclear.bam
	bash \$STAMPIPES/scripts/bam/random_sample.sh nuclear.bam subsample.bam 5000000
  samtools view -b -f 0x0040 subsample.bam > r1.bam

  # hotspot
  bash \$STAMPIPES/scripts/SPOT/runhotspot.bash \
    \$HOTSPOT_DIR \
    \$PWD \
    \$PWD/r1.bam \
    "${genome_name}" \
    "${params.readlength}" \
    DNaseI

  starch --header r1-both-passes/r1.hotspot.twopass.zscore.wig \
    > r1.spots.starch

  bash \$STAMPIPES/scripts/SPOT/info.sh \
    r1.spots.starch hotspot1 r1.spot.out \
    > r1.hotspot.info
  """
}

process bam_counts {
  label "modules"
  publishDir params.outdir

  input:
  path bam

  output:
  path 'tagcounts.txt'

  script:
  """
  python3 \$STAMPIPES/scripts/bwa/bamcounts.py \
    "${bam}" \
    tagcounts.txt
  """
}

process count_adapters {
  label "modules"
  publishDir params.outdir

  input:
  path bam

  output:
  path 'adapter.counts.txt'

  script:
  """
  bash "\$STAMPIPES/scripts/bam/count_adapters.sh" "${bam}" \
  | sed 's/^/adapter\t/' \
  > adapter.counts.txt
  """
}

process preseq {
  label "modules"
  publishDir params.outdir

  input:
  path nuclear_bam

  output:
  path 'preseq.txt'
  path 'preseq_targets.txt'
  path 'dups.hist'

  when:
  !params.UMI

  script:
  """
  python3 \$STAMPIPES/scripts/bam/mark_dups.py -i "${nuclear_bam}" -o /dev/null --hist dups.hist
  preseq lc_extrap -hist dups.hist -extrap 1.001e9 -s 1e6 -v > preseq.txt \
  || preseq lc_extrap -defects -hist dups.hist -extrap 1.001e9 -s 1e6 -v > preseq.txt

  # write out preseq targets
  bash "\$STAMPIPES/scripts/utility/preseq_targets.sh" preseq.txt preseq_targets.txt
  """
}

process cutcounts {
  label "modules"

  publishDir params.outdir

  label 'high_mem'

  input:
  path fai
  path filtered_bam

  output:
  path 'fragments.starch'
  path 'cutcounts.starch'
  path 'cutcounts.bw'
  path 'cutcounts.bed.bgz'
  path 'cutcounts.bed.bgz.tbi'

  script:
  """
  bam2bed --do-not-sort \
  < "${filtered_bam}" \
  | awk -v cutfile=cuts.bed -v fragmentfile=fragments.bed -f \$STAMPIPES/scripts/bwa/aggregate/basic/cutfragments.awk

  sort-bed fragments.bed | starch - > fragments.starch
  sort-bed cuts.bed | starch - > cuts.starch

  unstarch cuts.starch \
  | cut -f1-3 \
  | bedops -m - \
  | awk '{ for(i = \$2; i < \$3; i += 1) { print \$1"\t"i"\t"i + 1 }}' \
  > allbase.tmp

  unstarch cuts.starch \
  | bedmap --echo --count --delim "\t" allbase.tmp - \
  | awk '{print \$1"\t"\$2"\t"\$3"\tid-"NR"\t"\$4}' \
  | starch - > cutcounts.starch

  # Bigwig
  "\$STAMPIPES/scripts/bwa/starch_to_bigwig.bash" \
    cutcounts.starch \
    cutcounts.bw \
    "${fai}"

  # tabix
  unstarch cutcounts.starch | bgzip > cutcounts.bed.bgz
  tabix -p bed cutcounts.bed.bgz
  """
}

process density {
  label "modules"

  publishDir params.outdir
  label 'high_mem'

  input:
  path filtered_bam
  path chrom_bucket
  path fai

  output:
  path 'density.starch'
  path 'density.bw'
  path 'density.bgz'
  path 'density.bgz.tbi'
  tuple path(filtered_bam), path('density.starch'), emit: to_normalize

  script:
  def window_size = 75
  def bin_size = 20
  """
  mkfifo density.bed

  bam2bed -d \
  < "${filtered_bam}" \
  | cut -f1-6 \
  | awk '{ if( \$6=="+" ){ s=\$2; e=\$2+1 } else { s=\$3-1; e=\$3 } print \$1 "\t" s "\t" e "\tid\t" 1 }' \
  | sort-bed - \
  > density.bed \
  &

  unstarch "${chrom_bucket}" \
  | bedmap --faster --echo --count --delim "\t" - density.bed \
  | awk -v "binI=${bin_size}" -v "win=${window_size}" \
        'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print \$1 "\t" \$2 + shiftFactor "\t" \$3-shiftFactor "\tid\t" i \$4}' \
  | starch - \
  > density.starch

  # Bigwig
  "\$STAMPIPES/scripts/bwa/starch_to_bigwig.bash" \
    density.starch \
    density.bw \
    "${fai}" \
    "${bin_size}"

  # Tabix
  unstarch density.starch | bgzip > density.bgz
  tabix -p bed density.bgz
  """
}

process multimapping_density {

  publishDir params.outdir
  label 'modules'
  label 'high_mem'

  input:
  path marked_bam
  path chrom_bucket
  path fai

  output:
  path "mm_density.starch"
  path "mm_density.bw"
  path 'normalized.mm_density.starch'
  path 'normalized.mm_density.bw'

  script:
  def window_size = 75
  def bin_size = 20
  def scale = 1_000_000
  """
  # Mark multi-mapping reads as QC-pass!
  samtools view -h "${marked_bam}" |
  awk 'BEGIN{OFS="\t"} /XA:Z/ {\$2 = and(or(\$2, 2), compl(512))} 1' |
  samtools view --threads 3 -F 512 -o filtered.bam
  samtools index filtered.bam


  # Generate density
  mkfifo density.bed

  bam2bed -d \
  < filtered.bam \
  | cut -f1-6 \
  | awk '{ if( \$6=="+" ){ s=\$2; e=\$2+1 } else { s=\$3-1; e=\$3 } print \$1 "\t" s "\t" e "\tid\t" 1 }' \
  | sort-bed - \
  > density.bed \
  &

  unstarch "${chrom_bucket}" \
  | bedmap --faster --echo --count --delim "\t" - density.bed \
  | awk -v "binI=${bin_size}" -v "win=${window_size}" \
        'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print \$1 "\t" \$2 + shiftFactor "\t" \$3-shiftFactor "\tid\t" i \$4}' \
  | starch - \
  > mm_density.starch

  # Bigwig
  "\$STAMPIPES/scripts/bwa/starch_to_bigwig.bash" \
    mm_density.starch \
    mm_density.bw \
    "${fai}" \
    "${bin_size}"

  # # Tabix
  # unstarch density.starch | bgzip > density.bgz
  # tabix -p bed density.bgz

  rm density.bed

  # Normalized density
  unstarch mm_density.starch \
    | awk -v allcounts=\$(samtools view -c filtered.bam) \
          -v extranuclear_counts=\$(samtools view -c "filtered.bam" chrM chrC) \
          -v scale=${scale} \
          'BEGIN{ tagcount=allcounts-extranuclear_counts }
           { z=\$5;
             n=(z/tagcount)*scale;
             print \$1 "\t" \$2 "\t" \$3 "\t" \$4 "\t" n }' \
    | starch - > normalized.mm_density.starch

  "\$STAMPIPES/scripts/bwa/starch_to_bigwig.bash" \
    normalized.mm_density.starch \
    normalized.mm_density.bw \
    "${fai}" \
    "${bin_size}"
  """
}

process normalize_density {
  label "modules"
  publishDir params.outdir

  input:
  tuple path(filtered_bam), path(density)
  path fai

  output:
  path 'normalized.density.starch'
  path 'normalized.density.bw'
  path 'normalized.density.bgz'
  path 'normalized.density.bgz.tbi'

  script:
  def bin_size = 20
  def scale = 1_000_000
  """
  samtools index "${filtered_bam}"
  # Normalized density
  unstarch density.starch \
    | awk -v allcounts=\$(samtools view -c ${filtered_bam}) \
          -v extranuclear_counts=\$(samtools view -c "${filtered_bam}" chrM chrC) \
          -v scale=${scale} \
          'BEGIN{ tagcount=allcounts-extranuclear_counts }
           { z=\$5;
             n=(z/tagcount)*scale;
             print \$1 "\t" \$2 "\t" \$3 "\t" \$4 "\t" n }' \
    | starch - > normalized.density.starch

  "\$STAMPIPES/scripts/bwa/starch_to_bigwig.bash" \
    normalized.density.starch \
    normalized.density.bw \
    "${fai}" \
    "${bin_size}"

  unstarch normalized.density.starch | bgzip > normalized.density.bgz
  tabix -p bed normalized.density.bgz
  """
}

process insert_sizes {
  label "modules"

  publishDir params.outdir

  scratch false

  input:
  path nuclear_bam
  path nuclear_chroms

  output:
  path 'CollectInsertSizeMetrics.picard*'

  script:
  """
  export JAVA_TOOL_OPTIONS="-Djdk.lang.Process.launchMechanism=vfork"
  picard CollectInsertSizeMetrics \
    "INPUT=${nuclear_bam}" \
    OUTPUT=CollectInsertSizeMetrics.picard \
    HISTOGRAM_FILE=CollectInsertSizeMetrics.picard.pdf \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true

  cat CollectInsertSizeMetrics.picard \
  | awk '/## HISTOGRAM/{x=1;next}x' \
  | sed 1d \
  > hist.txt

  python3 "\$STAMPIPES/scripts/utility/picard_inserts_process.py" hist.txt > CollectInsertSizeMetrics.picard.info
  """
}

process motif_matrix {
  label "modules"

  publishDir params.outdir

  input:
  path hotspot_calls
  path fimo_transfac
  path fimo_names

  output:
  path 'hs_motifs*.txt'

  when:
  params.domotifs

  script:
  """
  # create sparse motifs
  bedmap --echo --echo-map-id --fraction-map 1 --delim '\t' "${hotspot_calls}" "${fimo_transfac}" > temp.bedmap.txt
  python \$STAMPIPES/scripts/bwa/aggregate/basic/sparse_motifs.py "${fimo_names}" temp.bedmap.txt
  """
}

process closest_features {
  label "modules"

  publishDir params.outdir

  input:
  path hotspot_calls
  path transcript_starts
  val thresholds

  output:
  path 'prox_dist.info'

  when:
  params.dofeatures

  script:
  """
  closest-features \
    --dist \
    --delim '\t' \
    --closest \
    "${hotspot_calls}" \
    "${transcript_starts}" \
  > closest.txt
  cat closest.txt \
  | grep -v "NA\$" \
  | awk -F"\t" '{print \$NF}' \
  | sed -e 's/-//g' \
  > closest.clean.txt

  for t in ${thresholds} ; do
    awk \
      -v t=\$t \
      '\$1 > t {sum+=1} END {print "percent-proximal-" t "bp " sum/NR}' \
      closest.clean.txt \
    >> prox_dist.info
  done
  """
}

/*
 * Metrics: Hotspot differential index comparison
 */
process differential_hotspots {
  label "modules"

  publishDir params.outdir

  input:
  path bam
  path peaks
  path index

  output:
  path 'differential_index_report.tsv'

  when:
  params.hotspot_index != "."

  script:
  def version = (new File(params.hotspot_index)).getAbsoluteFile().getParentFile().getName()
  def diffName = "dhsindex_${version}_differential_peaks"
  def diffPerName = "dhsindex_${version}_differential_peaks_percent"
  def conName = "dhsindex_${version}_constitutive_peaks"
  def conPerName = "dhsindex_${version}_constitutive_peaks_percent"

  """
  set -e -o pipefail
  statOverlap=\$(bedops -e 1 "${peaks}" "${index}" | wc -l)
  statNoOverlap=\$(bedops -n 1 "${peaks}" "${index}" | wc -l)
  total=\$(unstarch "${peaks}" | wc -l)
  statPercOverlap=\$(echo "scale=3; \$statOverlap * 100.0/\$total" | bc -q)
  statPercNoOverlap=\$(echo "scale=3; \$statNoOverlap * 100.0/\$total" | bc -q)

  {
    echo -e "${diffName}\t\$statNoOverlap"
    echo -e "${diffPerName}\t\$statPercNoOverlap"
    echo -e "${conName}\t\$statOverlap"
    echo -e "${conPerName}\t\$statPercOverlap"
  } > differential_index_report.tsv
  """
}


/*
 * Footprint calling
 */
process learn_dispersion {
  label "footprints"
  publishDir params.outdir
  memory '8 GB'
  cpus 8

  input:
  path ref
  path bam
  //path bai
  path spots
  path bias

  output:
  tuple path('dm.json'), path(bam), path("${bam}.bai"), emit: dispersion
  path 'dm.json', emit: to_plot

  when:
  params.bias != ""

  script:
  """
  samtools index "${bam}"

  # TODO: Use nuclear file
  unstarch ${spots} \
  | grep -v "_random" \
  | grep -v "chrUn" \
  | grep -v "chrM" \
  > intervals.bed

  ftd-learn-dispersion-model \
    --bm ${bias} \
    --half-win-width 5 \
    --processors 8 \
    ${bam} \
    ${ref} \
    intervals.bed \
  > dm.json
  """.stripIndent()
}

process make_intervals {

  label "footprints"

  input:
  path starch

  output:
  path 'chunk_*', emit: intervals

  script:
  """
  unstarch "${starch}" \
  | grep -v "_random" \
  | grep -v "chrUn" \
  | grep -v "chrM" \
  | split -l "${params.chunksize}" -a 4 -d - chunk_
  """.stripIndent()
}

process compute_deviation {
  label "footprints"
  memory '8 GB'
  cpus 4

  input:
  tuple path(interval), path(dispersion), path(bam), path(bai)
  path bias
  path ref

  output:
  path 'deviation.out', emit: deviations

  script:
  """
  ftd-compute-deviation \
  --bm "${bias}" \
  --half-win-width 5 \
  --smooth-half-win-width 50 \
  --smooth-clip 0.01 \
  --dm "${dispersion}" \
  --fdr-shuffle-n 50 \
  --processors 4 \
  "${bam}" \
  "${ref}" \
  "${interval}" \
  | sort --buffer-size=8G -k1,1 -k2,2n \
  > deviation.out
  """.stripIndent()
}

process merge_deviation {
  label "footprints"
  memory "32 GB"
  cpus 1

  input:
  path 'chunk_*'

  output:
  path 'interval.all.bedgraph', emit: merged_interval

  when:
  params.bias != ""

  script:
  """
  echo chunk_*
  sort -k1,1 -k2,2n -S 32G -m chunk_* > interval.all.bedgraph
  """.stripIndent()
}

process working_tracks {
  label "footprints"
  memory '32 GB'
  cpus 1
  publishDir params.outdir

  input:
  path merged_interval

  output:
  path 'interval.all.bedgraph', emit: bedgraph
  path 'interval.all.bedgraph.starch'
  path 'interval.all.bedgraph.gz'
  path 'interval.all.bedgraph.gz.tbi'

  script:
  """
  sort-bed "${merged_interval}" | starch - > interval.all.bedgraph.starch
  bgzip -c "${merged_interval}" > interval.all.bedgraph.gz
  tabix -0 -p bed interval.all.bedgraph.gz
  """.stripIndent()
}

process compute_footprints {
  label "footprints"
  memory '8 GB'
  cpus 1
  publishDir params.outdir

  input:
  tuple path(merged_interval), val(threshold)

  output:
  path "interval.all.fps.${threshold}.bed.gz"
  path "interval.all.fps.${threshold}.bed.gz.tbi"

  script:
  """
  output=interval.all.fps.${threshold}.bed
  cat "${merged_interval}" \
  | awk -v OFS="\t" -v thresh="${threshold}" '\$8 <= thresh { print \$1, \$2-3, \$3+3; }' \
	| sort-bed --max-mem 8G - \
	| bedops -m - \
	| awk -v OFS="\t" -v thresh="${threshold}" '{ \$4="."; \$5=thresh; print; }' \
  > \$output

  bgzip -c "\$output" > "\$output.gz"
  tabix -0 -p bed "\$output.gz"
  """.stripIndent()
}

process plot_footprints {

  label "footprints"
  publishDir params.outdir

  input:
  path model
  path plot

  output:
  path "dispersion.*pdf"

  script:
  """
  "./${plot}" "${model}"
  """
}

process cram {
  publishDir params.outdir
  cpus params.cramthreads / 2
  label 'samtools'

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

process starch_to_bigbed {
  publishDir "${params.outdir}"
  label "modules"

  input:
  path starch_in
  path chrom_sizes_bed

  output:
  path 'peaks/nuclear.peaks.fdr0.001.bb'

  script:
  outfile = starch_in.name.replace("starch", "bb")
  """
  cut -f1,3 "${chrom_sizes_bed}" > chrom_sizes
  mkdir -p peaks
  unstarch "${starch_in}" | cut -f1-4 > temp.bed
  bedToBigBed temp.bed chrom_sizes "peaks/${outfile}"
  rm temp.bed
  """
}

workflow {
  dataDir = "${baseDir}/../../../data"
  genome_name = file(params.genome).baseName

  // Create input channels
  bams = Channel.from(
      params.bams.tokenize(',')
    )
    .map {
      file(it)
    }
    .collect()

  // Main workflow execution
  merge_bams(bams)
  dups(merge_bams.out)
  filter_bam(dups.out.marked_bam)
  filter_nuclear(filter_bam.out.filtered_bam, file("${params.genome}.nuclear.txt"))

  // Conditional processes
  if (params.peakcaller == "macs2") {
    macs2(filter_nuclear.out.nuclear_bam)
  }

  if (params.peakcaller == "hotspot2") {
    hotspot2(
      params.hotspot_id,
      filter_nuclear.out.nuclear_bam,
      file(params.mappable),
      file(params.chrom_sizes),
      file(params.centers),
    )
  }

  // Analysis processes
  spot_score(
    filter_bam.out.filtered_bam,
    file("${dataDir}/annotations/${genome_name}.K${params.readlength}.mappable_only.bed"),
    file("${dataDir}/annotations/${genome_name}.chromInfo.bed"),
    genome_name,
  )

  bam_counts(dups.out.marked_bam)
  count_adapters(dups.out.marked_bam)

  if (!params.UMI) {
    preseq(filter_nuclear.out.nuclear_bam)
  }

  cutcounts(
    file("${params.genome}.fai"),
    filter_bam.out.filtered_bam,
  )

  density(
    filter_bam.out.filtered_bam,
    file(params.chrom_bucket),
    file("${params.genome}.fai"),
  )

  multimapping_density(
    dups.out.marked_bam,
    file(params.chrom_bucket),
    file("${params.genome}.fai"),
  )

  normalize_density(
    density.out.to_normalize,
    file("${params.genome}.fai"),
  )

  insert_sizes(
    filter_bam.out.filtered_bam,
    file("${params.genome}.nuclear.txt"),
  )

  // Conditional analysis processes
  if (params.peakcaller == "hotspot2") {
    if (params.domotifs) {
      motif_matrix(
        hotspot2.out.hotspot_calls,
        file("${dataDir}/motifs/${genome_name}.fimo.starch"),
        file("${dataDir}/motifs/${genome_name}.fimo.transfac.names.txt"),
      )
    }

    if (params.dofeatures) {
      closest_features(
        hotspot2.out.hotspot_calls,
        file("${dataDir}/features/${genome_name}.CombinedTxStarts.bed"),
        "0 1000 2500 5000 10000",
      )
    }

    if (params.hotspot_index != ".") {
      differential_hotspots(
        dups.out.marked_bam,
        hotspot2.out.onepercent_peaks,
        file(params.hotspot_index),
      )
    }

    // Footprint analysis
    if (params.bias != "") {
      learn_dispersion(
        file("${params.genome}.fa"),
        filter_bam.out.filtered_bam,
        hotspot2.out.hotspot_calls_for_bias,
        file(params.bias),
      )

      make_intervals(hotspot2.out.hotspot_calls)

      intervals_combined = make_intervals.out.intervals
        .flatten()
        .combine(learn_dispersion.out.dispersion)

      compute_deviation(
        intervals_combined,
        file(params.bias),
        file("${params.genome}.fa"),
      )

      merge_deviation(
        compute_deviation.out.deviations.collect()
      )

      working_tracks(
        merge_deviation.out.merged_interval
      )

      thresholds = Channel.from(0.2, 0.1, 0.05, 0.01, 0.001, 0.0001)
      intervals_with_thresholds = merge_deviation.out.merged_interval.combine(thresholds)

      compute_footprints(intervals_with_thresholds)

      plot_footprints(
        learn_dispersion.out.to_plot,
        file("${baseDir}/plot_footprints.py"),
      )
    }

    starch_to_bigbed(
      hotspot2.out.onepercent_peaks,
      file(params.chrom_sizes),
    )
  }

  // CRAM conversion
  bams_to_cram_combined = filter_bam.out.filtered_bam.mix(dups.out.bams_to_cram_marked)
  cram(
    bams_to_cram_combined,
    file("${params.genome}.fa"),
    file("${params.genome}.fai"),
  )
}
