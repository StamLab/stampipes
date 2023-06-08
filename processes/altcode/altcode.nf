nextflow.enable.dsl=2

// Workflows
workflow {

  STAR_solo(
    [
      tokenize_read_files(params.r1), tokenize_read_files(params.r2),
      params.barcodes_r1, params.barcodes_r2, params.barcodes_r3,
      [78, 8],
      [48, 8],
      [10, 8],
    ],
    file("${params.genome_dir}/*")
  )

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

  module 'STAR/2.7.9a'
  publishDir "test-out"

  input:
    tuple(
      path(r1), path(r2), // r1 and r2 may each receive multiple files
      path(r1_barcodes), path(r2_barcodes), path(r3_barcodes),
      val(r1_barcode_pos), val(r2_barcode_pos), val(r3_barcode_pos)
    )
    path("ref/*")

  output:
    path "output/Aligned.out.cram*", emit: cram
    path "output/Solo.out", emit: solo_analysis


  script:
    // TODO: How do we dynamically determine this?
    // barcode_positions = "0_10_0_17 0_48_0_55 0_78_0_85"
    bc1_position = pos_to_str(*r1_barcode_pos)
    bc2_position = pos_to_str(*r2_barcode_pos)
    bc3_position = pos_to_str(*r3_barcode_pos)
    umi_position = pos_to_str(0, 10)

    //bc1_position = pos_to_str(78, 8)
    //bc2_position = pos_to_str(48, 8)
    //bc3_position = pos_to_str(10, 8)
    //umi_position = pos_to_str(0, 10)

    // TODO: Determine from environment?
    bam_sort_RAM = 32_000_000_000

    r1_files = join_list_commas(r1)
    r2_files = join_list_commas(r2)

    num_threads = 10

    """
    mkdir -p output
    mkfifo output/Aligned.out.bam
    STAR \
    --genomeDir "ref" \
    --readFilesIn "${r1_files}" "${r2_files}" \
    --soloType CB_UMI_Complex \
    --soloCBposition "${bc3_position}" "${bc2_position}" "${bc1_position}" \
    --soloCBwhitelist "${r3_barcodes}" "${r2_barcodes}" "${r1_barcodes}" \
    --soloUMIposition "${umi_position}" \
    --soloCBmatchWLtype 1MM \
    --soloUMIdedup 1MM_All \
    --soloFeatures Gene GeneFull SJ \
    --runThreadN "${num_threads}" \
    --limitBAMsortRAM "${bam_sort_RAM}" \
    --outSAMtype BAM Unsorted \
    --outSAMattributes NH HI AS nM CR CY UR UY sM \
    --outBAMcompression -1 \
    - outBAMsortingThreadN "${num_threads}" \
    --readFilesCommand zcat \
    --outFileNamePrefix output/ \
    --limitOutSJcollapsed 5000000 &

    samtools sort \
      --reference  /net/seq/data2/projects/prime_seq/cell_ranger_ref/GRCh38-2022-Altius-gencode.v39-build/Homo_sapiens.GRCh38.dna.primary_assembly.fa.modified \
      -o output/Aligned.out.cram \
      --output-fmt-option "version=3.0,level=7" \
      --threads "${num_threads}" \
      --write-index \
      -T "tmpsort" \
      output/Aligned.out.bam &

    wait
    rm output/Aligned.out.bam
    """
}
