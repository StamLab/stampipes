nextflow.enable.dsl=2

params.outdir = "output"
params.metadata = ""

/// Workflows

/// Default workflow
/// Processes a single sample
workflow {

  def find_matrices_in_dir = { m, dir -> 
    def mtx_files = []
    dir.eachFileRecurse { f -> 
      // Find mtx.gz files, excluding the SJ folder
      if (f.name.endsWith(".mtx.gz") && !(f.parent.parent.name == "SJ")) {
        def parent_path = f.parent
        def relative_path = dir.parent.toUri().relativize( parent_path.toUri() ).toString()
        // Record the file, barcode file, feature file, and the relative path to store the output in
        mtx_files.push([m, file(params.metadata), f, file("${f.parent}/barcodes.tsv.gz"), file("${f.parent}/features.tsv.gz"), relative_path])
      }
    }
    return mtx_files
  }

  def meta = [:]
  def ref_files = file("${params.genome_dir}/*")

  STAR_solo(
    [
      meta,
      tokenize_read_files(params.r1), tokenize_read_files(params.r2),
      params.barcode_r1_list, params.barcode_r2_list, params.barcode_r3_list,
      [params.barcode_r1_pos, params.barcode_r1_len],
      [params.barcode_r2_pos, params.barcode_r2_len],
      [params.barcode_r3_pos, params.barcode_r3_len],
      [params.umi_pos, params.umi_len],
      ref_files,
      params.genome_fasta,
    ],
  )

  STAR_solo.out.solo_analysis
  | flatMap { find_matrices_in_dir(it[0], it[1]) }
  | convert_to_h5ad

  STAR_solo.out.solo_analysis
  | map { [it[0], it[1], file(params.metadata)] }
  | summarize_stats
}

// Helper functions
def tokenize_read_files(input) {
  if (input in String) {
    return input.tokenize(",")
  }
  return input
}

def join_list_commas(input) {
  // TODO: Do we need to handle different parameter passing styles differently?
  return input.join(",")
}

def pos_to_str(start, length) {
  return "0_${start}_0_${start+length-1}"
}

/// Processing
/// This process creates the Aligned.out.cram file and STARsolo analysis results
process STAR_solo {

  publishDir params.outdir, mode: "copy"
  cpus 30
  memory "80 GB"
  //scratch false

  input:
    tuple(
      val(meta),
      path(r1), path(r2), // r1 and r2 may each receive multiple files
      path(r1_barcodes), path(r2_barcodes), path(r3_barcodes),
      val(r1_barcode_pos), val(r2_barcode_pos), val(r3_barcode_pos),
      val(umi_barcode_pos),
      path("ref/*"),
      path(genome_fasta),
    )

  output:
    tuple(val(meta), path("Aligned.out.cram*"), emit: cram)
    tuple(val(meta), path("Solo.out"), emit: solo_analysis)


  script:
    // barcode_positions = "0_10_0_17 0_48_0_55 0_78_0_85"
    bc1_position = pos_to_str(*r1_barcode_pos)
    bc2_position = pos_to_str(*r2_barcode_pos)
    bc3_position = pos_to_str(*r3_barcode_pos)
    umi_position = pos_to_str(*umi_barcode_pos)

    // TODO: Determine from environment?
    bam_sort_RAM = 32_000_000_000

    r1_files = join_list_commas(r1)
    r2_files = join_list_commas(r2)

    num_threads = 30

    """
    set -o monitor
    mkfifo Aligned.out.bam

    (STAR \
      --soloCellReadStats Standard \
      --clip3pAdapterSeq AAAAAAAAAA \
      --clip3pAdapterMMp 0.1 \
      --genomeDir "ref" \
      --readFilesIn "${r1_files}" "${r2_files}" \
      --soloType CB_UMI_Complex \
      --soloCBposition "${bc3_position}" "${bc2_position}" "${bc1_position}" \
      --soloCBwhitelist "${r3_barcodes}" "${r2_barcodes}" "${r1_barcodes}" \
      --soloUMIposition "${umi_position}" \
      --soloCBmatchWLtype 1MM \
      --soloUMIdedup 1MM_All \
      --soloFeatures Gene GeneFull SJ GeneFull_Ex50pAS GeneFull_ExonOverIntron \
      --soloMultiMappers Unique PropUnique Uniform Rescue EM \
      --runThreadN "${num_threads}" \
      --limitBAMsortRAM "${bam_sort_RAM}" \
      --outSAMtype BAM Unsorted \
      --outSAMattributes NH HI AS nM CR CY UR UY sM \
      --outBAMcompression 0 \
      --outBAMsortingThreadN "${num_threads}" \
      --readFilesCommand zcat \
      --outFileNamePrefix ./ \
      --limitOutSJcollapsed 5000000 \
    || kill 0) &

    samtools sort \
      --reference  "${genome_fasta}" \
      -o Aligned.out.cram \
      --output-fmt-option "version=3.0,level=7" \
      --threads "${num_threads}" \
      --write-index \
      -T "tmpsort" \
      Aligned.out.bam

    wait
    rm Aligned.out.bam
    compress_mtx_files.sh ./Solo.out "${num_threads}"
    """
}

process convert_to_h5ad {
  cpus 1
  memory "10 GB"
  publishDir params.outdir, mode: "copy", saveAs: {f -> "$out_dir/$f"}

  input: 
    tuple(val(meta), path(metadata_file), path(matrix), path(barcodes), path(features), val(out_dir))

  output:
    tuple(val(meta), path(out_file))

  shell:
  out_file = "${matrix.simpleName}.h5ad"
  // scanpy requires specific file names
  '''
  mkdir -p tmp
  cp "!{matrix}" tmp/matrix.mtx.gz
  cp "!{barcodes}" tmp/barcodes.tsv.gz
  cp "!{features}" tmp/features.tsv.gz
  mtx_to_h5.py tmp "!{out_file}"  --metadata "!{metadata_file}" 
  rm -r tmp
  '''
}

process summarize_stats {
  cpus 1
  publishDir params.outdir, mode: "copy"

  input:
    tuple val(meta), path(solo_dir), path(metadata_json)

  output:
    tuple val(meta), path(stats_json)

  script:
    stats_json = "stats.json"
    """
    summarize_stats.py "${solo_dir}" "${metadata_json}" > "${stats_json}"
    """
}
