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
