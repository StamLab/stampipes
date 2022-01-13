nextflow.enable.dsl=2
params.outdir = "output"

params.id = ""

params.r1 = null
params.r2 = null

params.adapter_p5 = null
params.adapter_p7 = null
params.readlength = null
params.umimethod = null

params.starIndexDir = "/net/seq/data/genomes/human/GRCh38/noalts-sequins/STARgenome-gencode-v25"
ref_fasta = "${params.starIndexDir}/all.genome.fa.gz"
params.star_threads = 8


include { star } from "./modules/star.nf" addParams(publish: false)
include { fastp_adapter_trim } from "../../modules/adapter_trimming.nf"
include { move_umt; takara_trim_umt } from "../../modules/umt.nf"
include { publish } from "../../modules/utility.nf"
include { encode_cram; encode_cram_no_ref } from "../../modules/cram.nf" addParams(cram_write_index: false )

/// normalize_string_param coerces non-string types to the intended string type
/// and normalizes casing issues
def normalize_string_param(p) {
  switch (p) {
    case true:  // Regrettably, '--param ""' sets `param` to `true`.
    case false: // ?
    case null:  //
      return ""
    default:
        return p.toString().toLowerCase()
  }
}

workflow STAR_ALIGNMENT {
  
  main:
    def meta = [ id: params.id ]

    ref_files = file("${params.starIndexDir}/*")
    fastp_adapter_trim( [params.r1, params.r2, params.adapter_p5, params.adapter_p7] )

    // Decide which UMI filtering to use, if any
    switch (normalize_string_param(params.umimethod)) {
      case "takara-umt":
        takara_trim_umt(adapter_trim.out.fastq, params.readlength)
        star( takara_trim_umt.out.fastq, ref_files )
        move_umt(star.out.flatten())
        .set { aligned_bams }
        break

      case "":
      case "none":
        star( adapter_trim.out.fastq, ref_files )
        .set { aligned_bams }
        break

      default:
        throw new Exception("Can't understand umimethod parameter '${params.umimethod}'")
    }

    // Convert to CRAM
    // Transcriptome bams can't use the genome reference
    aligned_bams
      .flatten()
      .branch {
        transcriptome: it.name ==~ /.*toTranscriptome.*/
        genome: true
      }.set { bams }
    bams.genome.map {[meta, it, ref_fasta]} | encode_cram
    bams.transcriptome.map {[meta, it]} | encode_cram_no_ref

    // Publish output files
    encode_cram.out.cram
      .mix(encode_cram_no_ref.out.cram)
      .map { m, cram -> cram }
    | publish

}

workflow {
  STAR_ALIGNMENT()
}
