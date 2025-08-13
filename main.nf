// main.nf

// Define the pipeline schema and print help message
nextflow.enable.dsl=2

// Import modules and subworkflows
include { FASTQC } from './modules/nf-core/fastqc/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { ALIGNMENT_AND_QC } from './subworkflows/local/alignment_and_qc'
include { GATK_HAPLOTYPECALLER } from './modules/local/gatk_haplotypecaller'
include { VEP } from './modules/local/vep'

workflow {
    //
    // 1. Input channel and validation
    //
    Channel
        .fromPath(params.input)
        .splitCsv(header:true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.status = row.status
            [ meta, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ]
        }
        .set { ch_raw_reads }

    //
    // 2. Initial read quality control
    //
    FASTQC ( ch_raw_reads )

    //
    // 3. Alignment, sorting, and duplicate marking
    //
    ALIGNMENT_AND_QC ( ch_raw_reads )

    //
    // 4. Variant Calling
    //
    GATK_HAPLOTYPECALLER ( ALIGNMENT_AND_QC.out.bam_indexed )

    //
    // 5. Variant Annotation
    //
    VEP ( GATK_HAPLOTYPECALLER.out.vcf )

    //
    // 6. Aggregate results with MultiQC
    //
    // Create a channel with all the files for MultiQC
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT_AND_QC.out.samtools_stats)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT_AND_QC.out.markdup_metrics)
    ch_multiqc_files = ch_multiqc_files.mix(VEP.out.summary)

    MULTIQC (
        ch_multiqc_files.collect(),
        [] // empty meta map
    )
}
