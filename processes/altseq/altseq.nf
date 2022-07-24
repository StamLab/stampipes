nextflow.enable.dsl=2

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import com.xlson.groovycsv.CsvParser

include { BCL2DEMUX } from "../../modules/bcl2fastq.nf"
include { sort_and_encode_cram } from "../../modules/cram.nf"
include { publish_and_rename; publish } from "../../modules/utility.nf"

params.sample_config_tsv = ""
params.input_directory = ""


// Functions
def parse_sample_config(sample_config_tsv) {
  def data = new CsvParser().parseCsv(sample_config_tsv, separator: "\t")
  def sample_info = []
  for (sample in data) {
    sample_info.add( sample )
  } 
  return sample_info
}

def validate_sample_config(sample_info) {
  // Check for required info
  sample_info.collect {
    assert it.lane > 0 : "Sample has no lane: ${it}"
    assert it.barcode_index.size() > 0 : "Sample has no barcode index: ${it}"
    assert it.name : "Sample has name: ${it}"
  }
}


// test workflow
workflow test {

  def star_exe = file("${workflow.projectDir}/../../third_party/STAR")
  def genome_dir = file("/net/seq/data2/projects/prime_seq/cell_ranger_ref/star_2.7.10_genome_2022_gencode.v39/")
  def genome_fa = file("/net/seq/data2/projects/prime_seq/cell_ranger_ref/GRCh38-2022-Altius-gencode.v39-build/Homo_sapiens.GRCh38.dna.primary_assembly.fa.modified")
  def barcode_whitelist = file("/net/seq/data2/projects/prime_seq/barcodes-combined.txt")

  def sample_info = parse_sample_config(file(params.sample_config_tsv).text)
  validate_sample_config(sample_info)
  
  BCL2DEMUX(
    params.input_directory,
    sample_info,
    "s_1_1101"
  )
  | flatMap { it.sort(); it.collate(2) }
  | map { [
    sample_info.find(s -> it[0].baseName.startsWith("${s.name}_S") ),
    it[0],
    it[1],
  ]}
  | filter { it[0] != null }
  | view { "it is $it" }
  | set { fq_files }

  fq_files.flatMap { [ it[1], it[2] ] } | publish

  align(
    star_exe,
    genome_dir,
    barcode_whitelist,
    fq_files,
  )

  align.out.aligned_bam
  | map { [
    [name: it[0].name, id: it[0].name, barcode_index: it[0].barcode_index, lane: it[0].lane] ,
    it[1],
    genome_fa,
  ] }
  | sort_and_encode_cram

  sort_and_encode_cram.out.cram
  | map { ["${it[0].name}.sorted.cram", it[1]] }
  | publish_and_rename

}


// TODO once we figure out how it works
workflow ALTSEQ {

  }

process align {

  memory "108681M"
  cpus 5

  input:
    file star_exe
    file genome_dir
    file barcode_whitelist
    tuple val(meta), file(fq1), file(fq2)


  output:
    tuple val(meta), file("Aligned.out.bam"), emit: aligned_bam
    tuple val(meta), file("Solo.out"), emit: solo_directory

  shell:
    '''
    tmpdir=$(mktemp -d)
    "./!{star_exe}"   \
      --genomeDir "!{genome_dir}"   \
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
      --soloCellFilter  EmptyDrops_CR 96 .99 10 45000 90000 100000 0.01 20000 0.01 10000   \
      --quantMode "TranscriptomeSAM"   \
      --soloFeatures Gene GeneFull GeneFull_ExonOverIntron GeneFull_Ex50pAS   \
      --soloMultiMappers Unique PropUnique Uniform Rescue EM   \
      --readFilesCommand zcat   \
      --runThreadN 5   \
      --outSAMtype BAM Unsorted   \
      --outSAMattributes NH HI AS NM MD CR CY UR UY GX GN   \
      --outSAMunmapped Within   \
      --limitOutSJcollapsed 5000000   \
      --outTmpDir "$tmpdir/STARSolo"
    '''
}

