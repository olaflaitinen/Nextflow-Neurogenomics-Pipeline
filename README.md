# Nextflow-Neurogenomics-Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-239f20.svg)](https://www.nextflow.io/)
[![CI](https://github.com/your-username/Nextflow-Neurogenomics-Pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/your-username/Nextflow-Neurogenomics-Pipeline/actions/workflows/ci.yml)
[![nf-core](https://img.shields.io/badge/nf--core-compliant-brightgreen)](https://nf-co.re)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A `nf-core` compliant Nextflow pipeline for processing raw FASTQ files from whole-genome sequencing (WGS) of neurodegeneration patient cohorts. It performs read alignment, quality control, variant calling, and annotation.

### Pipeline Summary

The pipeline follows GATK best practices for germline short variant discovery.

```
                  ┌───┐
Input FASTQ files │...│
                  └───┘
                    │
                    ▼
           ┌────────────────┐
           │   FASTQC       │  (Read Quality Control)
           └────────────────┘
                    │
                    ▼
           ┌────────────────┐
           │    BWA-MEM     │  (Alignment to Reference Genome)
           └────────────────┘
                    │
                    ▼
           ┌────────────────┐
           │ SAMTOOLS SORT  │  (Sort & Index BAM)
           └────────────────┘
                    │
                    ▼
           ┌────────────────┐
           │ MARKDUPLICATES │  (GATK: Mark PCR Duplicates)
           └────────────────┘
                    │
                    ▼
           ┌────────────────┐
           │ HAPLOTYPECALLER│  (GATK: Variant Calling -> gVCF)
           └────────────────┘
                    │
                    ▼
           ┌────────────────┐
           │      VEP       │  (Variant Effect Predictor Annotation)
           └────────────────┘
                    │
                    ▼
           ┌────────────────┐
           │    MULTIQC     │  (Aggregate QC Report)
           └────────────────┘
```

### Quick Start

1.  **Install Nextflow (`>=21.10.3`)**

    ```bash
    curl -s https://get.nextflow.io | bash
    ```

2.  **Install Docker or Singularity**

    This pipeline requires a container engine. Docker is recommended for local execution.

3.  **Run the pipeline with test data**

    ```bash
    nextflow run your-username/Nextflow-Neurogenomics-Pipeline -r main -profile test,docker
    ```

### Usage

The pipeline requires a samplesheet and a reference genome.

**1. Samplesheet (`--input`)**

The input file must be a CSV file with the following columns: `sample`, `fastq_1`, `fastq_2`, `status`.

`assets/samplesheet.csv`:
```csv
sample,fastq_1,fastq_2,status
CONTROL01,./data/C01_R1.fastq.gz,./data/C01_R2.fastq.gz,control
PATIENT01,./data/P01_R1.fastq.gz,./data/P01_R2.fastq.gz,case
```

**2. Full Command**

```bash
nextflow run your-username/Nextflow-Neurogenomics-Pipeline \
    -profile docker \
    --input 'path/to/your/samplesheet.csv' \
    --genome 'GRCh38' \
    --outdir './results'
```

For more details on parameters and output, see `docs/USAGE.md` and `docs/OUTPUT.md`.

### Credits

This pipeline was built using the excellent `nf-core` community guidelines and tools. It uses the following software:

*   FastQC
*   BWA
*   Samtools
*   GATK4
*   VEP (Variant Effect Predictor)
*   MultiQC

Developed by O. Yunus L. Imanov.
