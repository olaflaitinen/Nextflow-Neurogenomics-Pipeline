# Nextflow-Neurogenomics-Pipeline: Usage

This document provides more detailed instructions on how to run the pipeline and explains the available parameters.

## Table of Contents

-   [Introduction](#introduction)
-   [Running the pipeline](#running-the-pipeline)
    -   [Quick Start](#quick-start)
    -   [Required Arguments](#required-arguments)
    -   [Reference Genomes](#reference-genomes)
    -   [Optional Arguments](#optional-arguments)
-   [Job Execution](#job-execution)

## Introduction

This pipeline automates the analysis of Whole-Genome Sequencing (WGS) data, from raw reads to annotated variants, following GATK best practices.

## Running the pipeline

### Quick Start

```bash
nextflow run olaflaitinen/Nextflow-Neurogenomics-Pipeline -r main -profile test,docker
```

### Required Arguments

-   `--input <samplesheet.csv>`

    Path to the input samplesheet. This must be a CSV file with the columns `sample,fastq_1,fastq_2,status`.

-   `--genome <genome_id>`

    Specifies the reference genome to use. The pipeline is configured to use `nf-core/iGenomes`. Common IDs include `GRCh37` and `GRCh38`. See the [iGenomes config](https://github.com/nf-core/configs/blob/master/conf/igenomes.config) for a full list.

### Reference Genomes

The pipeline uses pre-configured reference genomes from `nf-core/iGenomes`. When you specify `--genome GRCh38`, Nextflow will automatically download the required index files (BWA, FASTA, etc.) from an S3 bucket.

You can also provide your own reference files using parameters like `--bwa_index` and `--fasta`.

### Optional Arguments

-   `--outdir <directory>`

    The directory where results will be saved. Defaults to `./results`.

-   `--vep_cache_version <version>`

    The version of the VEP cache to use (e.g., `105`).

-   `--vep_species <species>`

    The species for VEP annotation (e.g., `homo_sapiens`).

## Job Execution

The pipeline can be run on different systems using Nextflow profiles.

-   `-profile docker`: Executes jobs in Docker containers. Recommended for local runs.
-   `-profile singularity`: Executes jobs in Singularity containers. Common on HPC systems.
-   `-profile slurm`: Submits jobs to a SLURM scheduler. Requires custom configuration for your cluster.

Example for running on a SLURM cluster:

```bash
nextflow run . -profile slurm,singularity \
    --input samples.csv \
    --genome GRCh38
```
```

---

### `modules/local/gatk_haplotypecaller.nf`

```groovy
// modules/local/gatk_haplotypecaller.nf

process GATK_HAPLOTYPECALLER {
    tag "$meta.id"
    publishDir "$params.outdir/variant_calling", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.g.vcf.gz"), emit: vcf
    tuple val(meta), path("*.g.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}"
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \\
        -R $fasta \\
        -I $bam \\
        -O ${prefix}.g.vcf.gz \\
        -ERC GVCF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk: \$(gatk --version 2>&1 | grep "The GATK")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.g.vcf.gz
    touch ${meta.id}.g.vcf.gz.tbi
    touch versions.yml
    """
}
