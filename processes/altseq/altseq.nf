nextflow.enable.dsl=2

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import com.xlson.groovycsv.CsvParser

include { BCL2DEMUX } from "../../modules/bcl2fastq.nf"
include { sort_and_encode_cram } from "../../modules/cram.nf"

params.sample_config_tsv = ""
params.input_directory = ""
params.star_exe = "${workflow.projectDir}/../../third_party/STAR"

params.outdir = "output"
params.publishmode = "link"

params.skip_alignment = false


// Functions
def parse_sample_config(sample_config_tsv) {
  def data = new CsvParser().parseCsv(sample_config_tsv, separator: "\t")
  def sample_info = []
  for (sample in data) {
    sample_info.add(
      sample.columns.collectEntries { c, v -> [c, sample[v]] }
    )
  } 
  return sample_info
}

def validate_sample_config(sample_info) {
  // Check for required info
  sample_info.collect {
    assert it.lane > 0 : "Sample has no lane: ${it}"
    assert it.barcode_index.size() > 0 : "Sample has no barcode index: ${it}"
    assert it.pool_name : "Sample has no pool name: ${it}"
    assert it.sample_name : "Sample has no sample name: ${it}"
  }
}


workflow ALTSEQ {

  // TODO: This list can probably be refined or repackaged!
  take:
    genome_dir
    genome_fa
    barcode_whitelist
    tiles
    input_dir
    sample_info

  main:

    // We demux at the library pool level.
    // So here, we do some work to get just those pools
    // We take the pool_name and just the first part of the barcode_index (the part before '-')
    sample_info_for_bcl2fastq = sample_info.collect { info -> [
      name: info.pool_name,
      lane: info.lane,
      barcode_index: info.barcode_index.split('-')[0],
    ]}.unique() // And filter to just the unique ones

    // Run BCL2Fastq on all files
    BCL2DEMUX(
      input_dir,
      sample_info_for_bcl2fastq,
      tiles,
    )

    // Merge and publish fastq files
    BCL2DEMUX.out
    | flatten()
    | map { fq_file ->
      // Extract the pool name, R1/R2, and the lane
      // output is shaped like [[prefix, read], [lane, file]]
      bcl_fq_regex = /(.*)_S[0-9]+_(L[0-9]+)_(R[1-2])_.*/
      match = (fq_file.baseName =~ bcl_fq_regex)[0]
      [ match[1,3], [match[2], fq_file]]
    }
    | filter {
      // We don't want to do further processing on Undetermined samples
      it[0][0] != "Undetermined"
    } 
    // Now we group it together by pool name, lane, and read 
    | groupTuple
    | map {
      readname, files -> [
        // Convert readname to string
        "${readname[0]}_${readname[1]}",
        //Make sure files are in order by lane
        files.sort { lane, filename -> lane <=> lane }
             .collect { lane, filename -> filename }
    ]}
    | merge_fq

    merge_fq.out
    // Use groupTuple to group files in R1, R2 pairs
    | map { [ it.baseName.replaceAll(/_R[12]/, "_RX"), it ] }
    | groupTuple(size: 2, sort: { a, b -> { a.baseName <=> b.baseName } } )
    // drop prefix, no longer needed
    | map { prefix, info -> info }
    // Re-associate the metadata
    | map {r1, r2 -> [
      sample_info_for_bcl2fastq.find(s -> r1.baseName.startsWith("${s.name}_R") ),
        r1,
        r2,
    ]}
    | set { merged_fq_files }

    if (!params.skip_alignment) {

      // Invoke STAR Solo
      align(
        genome_dir,
        genome_fa,
        params.star_exe,
        barcode_whitelist,
        merged_fq_files,
      )

      // Gather statistics and summarize the results

      // First, we pair up the analysis with the expected list of samples
      // (This key will help us decode pool/barcode -> sample)
      create_sample_configs(file(params.sample_config_tsv))
      | flatten()
      | map { fn -> [fn.baseName, fn] }
      | set { per_pool_sample_configs }

      align.out.solo_directory
      | map { meta, solo_dir -> ["${meta.name}_lane${meta.lane}", meta, solo_dir ] }
      | join(per_pool_sample_configs)
      | map { key, meta, solodir, config -> [meta, config, solodir] }
      | set {to_analyze}

      // Create pool-level stats and files
      analyze_solo_dir(to_analyze)

      // Create flowcell stats
      analyze_solo_dir.out.json_counts
      | map {meta, counts -> counts}
      | toSortedList
      | merge_stats

    }

    // Debugging section - use `nextflow run -dump-channels` to write channel contents to terminal
    //merged_fq_files.dump(tag: "merged_fq_files")
    //per_pool_sample_configs.dump(tag: "per_pool")
    //to_analyze.dump(tag: "to_analyze")
}

workflow {

  def genome_dir = file(params.genome_dir)
  def genome_fa = file(params.genome_fa)
  def barcode_whitelist = file(params.barcode_whitelist)

  def sample_info = parse_sample_config(file(params.sample_config_tsv).text)
  validate_sample_config(sample_info)

  ALTSEQ(genome_dir, genome_fa, barcode_whitelist, "s_*", params.input_directory, sample_info)
}

// test workflow
 workflow test {
  println "Running test workflow..."

  def star_exe = file("${workflow.projectDir}/../../third_party/STAR")
  def genome_dir = file("/net/seq/data2/projects/prime_seq/cell_ranger_ref/star_2.7.10_genome_2022_gencode.v39/")
  def genome_fa = file("/net/seq/data2/projects/prime_seq/cell_ranger_ref/GRCh38-2022-Altius-gencode.v39-build/Homo_sapiens.GRCh38.dna.primary_assembly.fa.modified")
  def barcode_whitelist = file("/net/seq/data2/projects/prime_seq/barcodes-combined.txt")

  def sample_info = parse_sample_config(file(params.sample_config_tsv).text)
  validate_sample_config(sample_info)

  ALTSEQ(genome_dir, genome_fa, barcode_whitelist, "s_[1-4]_1234", params.input_directory, sample_info)
}


process align {

  memory "108681M"
  cpus {cpus}
  scratch false  // Was filling up tmp dirs
  tag "${meta.name}"

  input:
    path genome_dir
    path reference_fa
    path star_exe
    path barcode_whitelist
    tuple val(meta), path(fq1), path(fq2)


  output:
    tuple val(meta), file("Aligned.out.cram"), file("Aligned.out.cram.crai"), emit: aligned_cram
    tuple val(meta), file("Solo.out"), emit: solo_directory
    tuple val(meta), file("Log.out"), file("Log.final.out"), emit: logs

  shell:
    cpus = 5
    cram_fmt_options = [
      "version=3.0",
      "level=7",
      "lossy_names=0",
      ].join(",")
    bam_fifo = "Aligned.out.bam"
    '''
    set -e
    tmpdir=$(mktemp -d)
    echo "Tmpdir is: ${tmpdir}"
    mkfifo "!{bam_fifo}"
    "./!{star_exe}"   \
      --genomeDir "!{genome_dir}"   \
      --genomeLoad LoadAndRemove \
      --readFilesIn "!{fq2}" "!{fq1}"   \
      --soloType CB_UMI_Simple   \
      --soloCellReadStats Standard   \
      --clip3pAdapterSeq AAAAAAAAAA   \
      --clip3pAdapterMMp 0.1   \
      --soloCBstart 1   \
      --soloCBlen 12   \
      --soloUMIstart 13   \
      --soloUMIlen 16   \
      --soloCBwhitelist "!{barcode_whitelist}"   \
      --soloFeatures GeneFull_Ex50pAS   \
      --soloMultiMappers Unique PropUnique Uniform Rescue EM \
      --readFilesCommand zcat   \
      --runThreadN "!{cpus}"   \
      --outSAMtype BAM Unsorted   \
      --outBAMcompression 0 \
      --outSAMattributes NH HI AS NM MD CR CY UR UY GX GN \
      --outSAMunmapped Within   \
      --limitOutSJcollapsed 5000000   \
      --outTmpDir "$tmpdir/STARSolo" \
      &  # Launch in background, so we can convert to cram from pipe

    samtools sort \
      --reference "!{reference_fa}" \
      --output-fmt-option "!{cram_fmt_options}" \
      --threads "!{cpus}" \
      -@ "!{cpus}" \
      -o Aligned.out.cram \
      --write-index \
      -m 2G \
      -M \
      "!{bam_fifo}"


    wait # probably not necessary
    rm -rf "$tmpdir" # FIXME: This doesn't get called if STAR crashes or if NF cancels
    if [[ $(wc -l Solo.out/GeneFull_Ex50pAS/CellReads.stats) -le 2 ]] ; then
      echo -e "CellReads.stats does not contain enough output.\n" \
              "This may be caused by running out of /tmp space.\n" \
              "Alternately, it could be that no matching barcodes were found.\n" \
              >&2 
      exit 1
    fi
    '''
}


process merge_fq {

  cpus {cpus}
  container null
  module "htslib/1.12"
  scratch false

  input:
    tuple val(name), path("in.*.fq.gz")

  output:
    file out

  shell:
    cpus = 10
    out = "${name}.fq.gz"
    '''
    zcat in.*.fq.gz \
    | bgzip --stdout --threads "!{cpus}" \
    > "!{out}"
    '''
}

process analyze_solo_dir {
  scratch false

  input:
    tuple val(meta), path(sample_config), path("Solo.out")

  output:
    tuple val(meta), file("output")
    tuple val(meta), file("output/counts.json"), emit: "json_counts"

  shell:
    '''
    sed 's/[ACTGN]*-//' < '!{sample_config}' > barcode.config
    for dir in GeneFull_Ex50pAS ; do
      outdir=output/$dir
      allcountsfile=$outdir/allcounts.csv
      mkdir -p "$outdir"
      bash matrix2csv.sh  "Solo.out/$dir/raw/" > "$allcountsfile"
      cat barcode.config | while read name barcode ; do
        cat "$allcountsfile" \
        | awk -F, -vbarcode=$barcode -vname=$name \
          '$1 == barcode { print $2 "," $3 "," $5 }' \
          > "$outdir/$name.counts.csv"
      done
      analyze.py "Solo.out/$dir/CellReads.stats" "barcode.config" "$outdir"
    done
    generate_counts_json.py "Solo.out/GeneFull_Ex50pAS" "!{sample_config}" "!{meta.name}" > output/counts.json
    '''
}

process create_sample_configs {
  scratch false
  executor "local"
  // TODO: Take sample_config as val
  // That way we don't rely on column ordering
  input:
    path sample_config
  output:
    file("configs/*")

  shell:
    '''
    mkdir configs
    awk < '!{sample_config}' \
    'NR > 1 { print $2 "\t" $4 > "configs/" $1 "_lane" $3 ".config"}'
    '''
}

process merge_stats {
  scratch false
  executor "local"
  input:
    path("*.input.json")
  output:
    path("flowcell_stats.json")
  shell:
  '''
  jq --slurp --sort-keys -c . *.input.json > flowcell_stats.json
  '''
}
