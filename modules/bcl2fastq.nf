nextflow.enable.dsl=2

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import com.xlson.groovycsv.CsvParser

params.input_directory = ""
params.bcl2fastq_threads = 20

params.bcl2fastq_sample_config_tsv = ""
params.bases_mask = ""

params.bcl2fastq_tiles = "*"


params.samplesheet_header = """[Header]
                              |Workflow,GenerateFASTQ""".stripMargin()

// Functions

def parse_sample_config(sample_config_tsv) {
  def data = new CsvParser().parseCsv(sample_config_tsv, separator: "\t")
  def sample_info = []
  for (sample in data) {
    sample_info.add( sample )
  } 
  //sample_info = sample_info.collect { it[]}
  for (s in sample_info) {
    println "sample is ${s}"
  }
  return sample_info
}

def validate_sample_config(sample_info) {
  // Check for required info
  sample_info.collect {
    // println "Sample info: ${it}"
    assert it.lane > 0 : "Sample has no lane: ${it}"
    assert it.barcode_index.size() > 0 : "Sample has no barcode index: ${it}"
    assert it.name : "Sample has name: ${it}"
  }
}


workflow {
  def txt = file(params.bcl2fastq_sample_config_tsv).text
  def sample_info = parse_sample_config(txt)
  validate_sample_config(sample_info)

  BCL2DEMUX(
    params.input_directory,
    Channel.from(sample_info),
    "*",
  )

}

workflow test {
  // Dunno yet how to test this... Maybe with a real flowcell, but only a couple of tiles of data?
  def tiles = "s_1_0002"
  def input_dir = "/net/seq/data2/sequencers/220627_A01698_0053_BH2TCJDSX5"
  def txt = """name\tbarcode_index\tlane
              |N701\tTAAGGCGA\t1
              |N702\tCGTACTAG\t1
              |N703\tAGGCAGAA\t1
              |N704\tTCCTGAGC\t1
              |N705\tGGACTCCT\t1
              |N706\tTAGGCATG\t1
              |N707\tCTCTCTAC\t1
              |N708\tCAGAGAGG\t1
              |N709\tGCTACGCT\t1
              |N710\tCGAGGCTG\t1
              |N711\tAAGAGGCA\t1
              |N712\tGTAGAGGA\t1""".stripMargin()

  def sample_info = parse_sample_config(txt)
  validate_sample_config(sample_info)


  BCL2DEMUX (
    input_dir,
    sample_info,
    tiles
  )

}


workflow BCL2DEMUX {

  take:
    illumina_dir 
    sample_info 
    tiles
    // TODO

  main:
    text = generate_samplesheet( [params.samplesheet_header, sample_info])
    

  //emit:
    // TODO
  
}


// process demux_python {

// }

process generate_samplesheet {

  input:
    tuple val(header), val(sample_info)

  output:
    file("Samplesheet.csv")

  shell:
    sheet = [
      header,
      "",
      "[Settings]",
      "Lane,SampleID,index",
      *sample_info.collect { "${it.lane},${it.name},${it.barcode_index}" }
    ].join("\n")


    '''
    printf '!{sheet}' > Samplesheet.csv
    '''
}

// process generate_samplesheet_nodemux {

// }


process bcl2fastq {

  container "dceoy/bcl2fastq@sha256:6d7233f2160721d6cb62f77a127d499597f4b35bb435cc8265d05f5bf54c7b94"

  input:
    tuple file(illumina_dir), file(samplesheet), val(tiles)

  output:
    file("output/*")

  shell:
  '''
    mkdir output
    bcl2fastq \
      --input-dir "!{illumina_dir}/Data/Intensities/BaseCalls" \
      --samplesheet "${samplesheet}" \
      --use-bases-mask "!{bcl_mask}" \
      --output-dir "output/" \
      --barcode-mismatches "!{mismatches}" \
      --tiles "!{tiles}" \
      --loading-threads        $(( SLURM_CPUS_PER_TASK / 2 )) \
      --writing-threads        $(( SLURM_CPUS_PER_TASK / 2 )) \
      --processing-threads     $(( SLURM_CPUS_PER_TASK ))
  '''
}

// process bclconvert {

// }
