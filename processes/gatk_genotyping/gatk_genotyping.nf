/// run bwa and GATK on trimmed reads


params.saveMode = 'copy'
// these variables are provided by the sample info or the processing bash script.
//params.flowcellID = ""
//params.sampleName = ""
//params.libraryID = ""
params.trimmomatic = "/net/seq/data2/projects/gchalla/WGS_Pilot/Trimmomatic-0.39"
params.genomeRefDir = "/net/seq/data2/projects/gchalla/Broad_Genome_Ref_v0"
params.genomeRef = "Homo_sapiens_assembly38.fasta"
params.gatk_exe = "/net/seq/data2/projects/gchalla/WGS_Pilot/gatk-4.3.0.0/gatk"
params.suffix = 'hs38'

params.inputDir = ""
params.outDir = ""
params.resultsDir = 'results/trimmomatic'
params.memResultsDir = 'results/bwa_mem'

params.haplotypeCallerResultsDir = 'results/gatk/haplotypeCaller'

params.filePattern = "./*_{R1,R2}.fastq.gz"
params.trimfilePattern = "./*_paired_{R1,R2}.fastq.gz"

Channel.fromFilePairs(params.filePattern)
        .into { ch_in_trimmomatic }

Channel.trimfromFilePairs(params.filePattern)
        .into { ch_in_mem }

Channel.value("$params.genomeRefDir/$params.rgenomeRef")
        .set { ch_refFasta }

workflow {
    TRIMMOMATIC
    bwamem
    bamsort
    markdup
    bqsr
    apply_bqsr
    haplotype_caller
    genotyping
}

process TRIMMOMATIC {

    publishDir params.resultsDir, mode: params.saveMode

    input:
    tuple flowcellID, libraryID, sampleName, file(reads) from ch_in_trimmomatic
    path refFasta from ch_refFasta

    output:
    tuple path(fq_1_paired), path(fq_2_paired) into ch_out_trimmomatic

    script:

    fq_1_paired = flowcellID + sampleName + libraryID + '_paired_R1.fastq.gz'
    fq_1_unpaired = flowcellID + sampleName + libraryID +'_singleton_R1.fastq.gz'
    fq_2_paired = flowcellID + sampleName + libraryID + '_paired_R2.fastq.gz'
    fq_2_unpaired = flowcellID + sampleName + libraryID + '_singleton_R2.fastq.gz'



    """
    module load jdk/11.0.16
    java -jar params.trimmomatic/trimmomatic-0.39.jar \
    PE -phred33 \
    -threads 8
    ${reads[0]} \
    ${reads[1]} \
    $fq_1_paired \
    $fq_1_unpaired \
    $fq_2_paired \
    $fq_2_unpaired \
    ILLUMINACLIP:${ADAPTERS}/NexteraPE-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """

}

process bwamem {

    publishDir params.memResultsDir, mode: params.saveMode

    input:
    tuple file(trimmedReads) from ch_out_trimmomatic
    set flowcellID, libraryID, sampleName, suffix

    output:
    file('*.bam') into ch_out_bwamem
    
    script:
    READGROUP="@RG\\tID:${flowcellID}_${sampleName}_${libraryID}\\tSM:$sampleName\\tLB:$libraryID\\tPL:ILLUMINA"

    """
    module load bwa/0.7.17 
    module load samtools/1.7
    cp ${params.indexResultsDir}/* .
    cp ${params.samtoolsFaidxResultsDir}/* .
    bwa mem -K 100000000 -Y  -R "${TAG}\" ${params.refFasta} ${trimmedReads[0]} ${trimmedReads[1]} \
        | samtools view -@ 8 -Sb -t ${params.refFasta} \
        -o ${flowcellID}_${sampleName}_${libraryID}_${suffix}_unsorted.bam
    """
}

process bamsort {
    publishDir params.bamsortResultsDir, mode: params.saveMode

    input:
    file('*.bam') from ch_out_bwamem

    output:
    file('*.bam') into ch_out_bamsort
    
    script:
    
    """
    module load samtools/1.7
    samtools sort -@ 8 -m 4G ${flowcellID}_${sampleName}_${libraryID}_${suffix}_unsorted.bam \
        -O BAM -o ${flowcellID}_${sampleName}_${libraryID}_${suffix}_sorted.bam
    """
}

process markdup {
    publishDir params.memResultsDir, mode: params.saveMode

    input:
    file('*.bam') from ch_out_bamsort

    output:
    file('*.bam') into ch_out_markdup
    
    script:
    
    """
        ${params.gatk} MarkDuplicates \
            -I ${flowcellID}_${sampleName}_${libraryID}_${suffix}_sorted.bam \
            -O ${flowcellID}_${sampleName}_${libraryID}_${suffix}_markdup.bam \
            --METRICS_FILE ${flowcellID}_${sampleName}_${libraryID}_${suffix}_markdup.log \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --CREATE_INDEX \
            -R ${params.refFasta}
    """
}

process bqsr {
    publishDir params.memResultsDir, mode: params.saveMode

    input:
    file('*.bam') from ch_out_markdup


    output:
    file('*_markdup_recal.table') into ch_out_bqsr
    
    script:
    
    """
        ${params.gatk} BaseRecalibrator \
            -I ${flowcellID}_${sampleName}_${libraryID}_${suffix}_markdup.bam \
            -O ${flowcellID}_${sampleName}_${libraryID}_${suffix}_markdup_recal.table \
            -R ${params.refFasta} \
            --known-sites ${params.genomeRefDir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            --known-sites ${params.genomeRefDir}/1000G_phase3_v4_20130502.sites.hg38.vcf \
            --known-sites ${params.genomeRefDir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz
    """
}
process apply_bqsr {
    publishDir params.memResultsDir, mode: params.saveMode

    input:
    file('*_markdup_recal.table') from ch_out_bqsr

    output:
    file('*.bam') into ch_out_applyBQSR
    
    script:
    
    """
        ${params.gatk} ApplyBQSR \
            -I ${flowcellID}_${sampleName}_${libraryID}_${suffix}_markdup.bam \
            -O ${flowcellID}_${sampleName}_${libraryID}_${suffix}_markdup_recal.bam \
            --bqsr-recal-file ${flowcellID}_${sampleName}_${libraryID}_${suffix}_markdup_recal.table \
            -R ${params.refFasta} \
            --create-output-bam-index
    """
}

process haplotype_caller {
    publishDir params.memResultsDir, mode: params.saveMode

    input:
    file('*.bam') from ch_out_applyBQSR

    output:
    file('*_raw_germline_haplocall.vcf.gz') into ch_out_haplocall
    
    script:
    
    """
        ${params.gatk} HaplotypeCaller \
            -I ${flowcellID}_${sampleName}_${libraryID}_${suffix}_markdup_recal.bam \
            -O ${flowcellID}_${sampleName}_${libraryID}_${suffix}_raw_germline_haplocall.vcf.gz \
            -ERC GVCF \
            -R ${params.refFasta} 
    """
}

process genotyping {
    publishDir params.memResultsDir, mode: params.saveMode

    input:
    file('*_raw_germline_haplocall.vcf.gz') from ch_out_haplocall

    output:
    file('*_raw_germline.vcf.gz') into ch_out_genotype
    
    script:
    
    """
        ${params.gatk} GenotypeGVCFs \
            -R ${params.refFasta} \
            -V ${flowcellID}_${sampleName}_${libraryID}_${suffix}_raw_germline_haplocall.vcf.gz \
            -O ${flowcellID}_${sampleName}_${libraryID}_${suffix}_raw_germline.vcf.gz
    """
}