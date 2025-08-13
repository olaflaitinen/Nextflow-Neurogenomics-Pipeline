// modules/local/samtools_sort.nf

process SAMTOOLS_SORT {
    tag "$meta.id"
    publishDir "$params.outdir/alignment", mode: 'copy', pattern: "*.{bam,bai}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    tuple val(meta), path("*.sorted.bam.bai"), emit: bai
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}.sorted.bam"
    """
    samtools sort -@ ${task.cpus - 1} -o $prefix $bam
    samtools index $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    touch ${meta.id}.sorted.bam.bai
    touch versions.yml
    """
}
