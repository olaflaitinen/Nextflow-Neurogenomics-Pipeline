// subworkflows/local/alignment_and_qc.nf

include { BWA_MEM } from '../../modules/local/bwa_mem'
include { SAMTOOLS_SORT } from '../../modules/local/samtools_sort'
// ... import other modules like MarkDuplicates, SamtoolsStats, etc.

workflow ALIGNMENT_AND_QC {
    take:
    ch_reads // channel: [ meta, [ fq1, fq2 ] ]

    main:
    ch_versions = Channel.empty()

    // Get reference files from iGenomes
    // (This logic would typically be in the main workflow or a conf file)
    ch_bwa_index = file("${params.genomes_base}/${params.genome}/BWAIndex/**")

    BWA_MEM ( ch_reads, ch_bwa_index )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    SAMTOOLS_SORT ( BWA_MEM.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    // ... chain other processes like MarkDuplicates here ...
    
    emit:
    bam_indexed = SAMTOOLS_SORT.out.bam_indexed // Final output of this subworkflow
    // samtools_stats = SAMTOOLS_STATS.out.stats
    // markdup_metrics = MARKDUP.out.metrics
    versions = ch_versions.ifEmpty(null)
}
