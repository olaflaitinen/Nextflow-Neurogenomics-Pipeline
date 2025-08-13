// modules/local/vep.nf

process VEP {
    tag "$meta.id"
    publishDir "$params.outdir/annotation", mode: 'copy'

    input:
    tuple val(meta), path(vcf)
    path vep_cache
    path fasta

    output:
    tuple val(meta), path("*.vep.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vep.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.html"), emit: summary
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}.vep"
    """
    vep --input_file $vcf \\
        --output_file ${prefix}.vcf \\
        --stats_file ${prefix}.summary.html \\
        --cache --dir_cache $vep_cache \\
        --fasta $fasta \\
        --format vcf --vcf \\
        --symbol --terms SO --tsl \\
        --hgvs --fasta $fasta \\
        --fork ${task.cpus} \\
        --offline

    bgzip ${prefix}.vcf
    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(vep --help | head -n 2 | tail -n 1 | sed 's/^#.*ensembl-vep/ensembl-vep/')
    END_VERSIONS
    """
}
