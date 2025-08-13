// modules/local/bwa_mem.nf

process BWA_MEM {
    tag "$meta.id"
    publishDir "$params.outdir/alignment", mode: 'copy', pattern: "*.bam"

    input:
    tuple val(meta), path(reads)
    path bwa_index

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    def fastq_1 = reads[0]
    def fastq_2 = reads[1]
    def prefix = "${meta.id}"

    // Use GATK-compatible read group header
    def read_group = "\'@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:WGS\'"

    """
    bwa mem -K 100000000 -p -v 3 -t $task.cpus -R $read_group \\
        ${bwa_index[0].baseName} \\
        $fastq_1 $fastq_2 | samtools view -bS - > ${prefix}.bam
    """
    
    stub:
    """
    touch ${meta.id}.bam
    """
}
