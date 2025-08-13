# Nextflow-Neurogenomics-Pipeline: Output

This document describes the output directory structure from the pipeline.

## Pipeline overview

The pipeline is built using Nextflow and processes data as follows:

-   Raw read QC (`FastQC`)
-   Alignment to reference (`BWA-MEM`)
-   Sorting and duplicate marking (`Samtools`, `GATK MarkDuplicates`)
-   Variant calling (`GATK HaplotypeCaller`)
-   Variant annotation (`Ensembl VEP`)
-   Aggregated QC (`MultiQC`)

## Output Directory Structure

The results will be created in the directory specified by `--outdir`.

<outdir>/
├── alignment/ # Alignment files
│ ├── <sample>.bam
│ └── <sample>.bam.bai
├── variant_calling/ # Called variants
│ ├── <sample>.g.vcf.gz
│ └── <sample>.g.vcf.gz.tbi
├── annotation/ # Annotated variants
│ ├── <sample>.annotated.vcf.gz
│ ├── <sample>.annotated.vcf.gz.tbi
│ └── <sample>.vep.summary.html
├── qc/ # Quality control files
│ ├── fastqc/ # Raw read QC
│ └── samtools_stats/ # Alignment QC
├── multiqc/ # Aggregated MultiQC report
│ ├── multiqc_report.html
│ └── multiqc_data/
└── pipeline_info/ # Nextflow execution reports
├── execution_report.html
└── ...
<outdir>

### Key Output Files

-   **`<outdir>/alignment/<sample>.bam`**: The final, sorted, duplicate-marked alignment file for each sample. This is the primary alignment output.
-   **`<outdir>/variant_calling/<sample>.g.vcf.gz`**: Per-sample variants in gVCF format, produced by GATK HaplotypeCaller.
-   **`<outdir>/annotation/<sample>.annotated.vcf.gz`**: The final VCF file with annotations added by VEP.
-   **`<outdir>/multiqc/multiqc_report.html`**: The main summary report that aggregates QC metrics from all steps and all samples. **This is the best place to start when reviewing results.**
