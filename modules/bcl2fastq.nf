nextflow.enable.dsl=2

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import com.xlson.groovycsv.CsvParser

params.input_directory = ""
params.bcl2fastq_threads = 20

params.bcl2fastq_sample_config_tsv = ""
params.bases_mask = ""

params.bcl2fastq_tiles = "s_*"


params.samplesheet_header = """[Header]
                              |Project Name,Stamlab
                              |Workflow,GenerateFASTQ""".stripMargin()

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


workflow {
  def txt = file(params.bcl2fastq_sample_config_tsv).text
  def sample_info = parse_sample_config(txt)
  validate_sample_config(sample_info)

  BCL2DEMUX(
    params.input_directory,
    Channel.from(sample_info),
    params.bcl2fastq_tiles,
  )

}

workflow test {
  // Test with a subset of tiles
  def tiles = "s_1_121[0-9]"
  def input_dir = "/net/seq/data2/sequencers/220627_A01698_0053_BH2TCJDSX5"
  def config_txt = """name\tbarcode_index\tlane
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

  def sample_info = parse_sample_config(config_txt)
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

  main:

    // Group samples by lanes of processing
    // This will let us run 1 bcl2fastq job for each lane
    Channel.fromList(sample_info)
    | map { [it.lane, it] }
    | groupTuple(sort: 'hash')
    | map { it[1] }
    | set { sample_info_by_lane }

    sample_info_by_lane
    | map { [params.samplesheet_header, it] }
    | generate_samplesheet
    | map { it -> [ illumina_dir, it, tiles] }
    | bcl2fastq

  emit:
    bcl2fastq.out
  
}


// process demux_python {

// }

process generate_samplesheet {
  executor = "local"

  input:
    tuple val(header), val(sample_info)

  output:
    file("Samplesheet.csv")

  // TODO: This is silly. We create a bash script to just
  // write a file.  However, despite my best efforts, I
  // haven't figured out how to get this to work as an
  // `exec` block.
  shell:
    settings = ""
    sheet_parts = [
      header,
      "",
      "[Settings]",
      settings,
      "[Data]",
      "Lane,SampleID,SampleName,index",
      *sample_info.collect { "${it.lane},${it.name},${it.name},${it.barcode_index}" }
    ]
    sheet = sheet_parts.join("\n")

    '''
    echo '!{sheet}' > Samplesheet.csv
    '''
}

// process generate_samplesheet_nodemux {

// }


process bcl2fastq {

  container "dceoy/bcl2fastq@sha256:6d7233f2160721d6cb62f77a127d499597f4b35bb435cc8265d05f5bf54c7b94"
  cpus {cpus}
  scratch false

  input:
    tuple path(illumina_dir), path(samplesheet), val(tiles)

  output:
    file("output/*fastq.gz")

  shell:
    cpus = 10
    '''
    workdir=$PWD
    outdir=$workdir/output
    mkdir -p "$outdir"
    cd "!{illumina_dir}"
    bcl2fastq \
      --input-dir "Data/Intensities/BaseCalls" \
      --sample-sheet "$workdir/!{samplesheet}" \
      --output-dir "$outdir" \
      --tiles "!{tiles}" \
      --barcode-mismatches 1
    '''
}

// process bclconvert {

// }
