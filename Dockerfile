# Start from a base image with common utilities
FROM ubuntu:20.04

# Set environment variables to avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV LANG C.UTF-8

# Install base dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    openjdk-11-jre \
    wget \
    curl \
    unzip \
    git \
    bzip2 \
    ca-certificates \
    perl \
    python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Install Samtools and BWA
RUN apt-get update && apt-get install -y samtools bwa

# Install FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    ln -s /FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.11.9.zip

# Install GATK4
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip && \
    unzip gatk-4.2.0.0.zip && \
    rm gatk-4.2.0.0.zip
# Add GATK to path
ENV PATH="/gatk-4.2.0.0:${PATH}"
# Create a wrapper script for GATK
RUN echo '#!/bin/bash\njava -jar /gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar "$@"' > /usr/local/bin/gatk && \
    chmod +x /usr/local/bin/gatk

# Install VEP
RUN wget https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/105.0.zip && \
    unzip 105.0.zip && \
    rm 105.0.zip && \
    cd ensembl-vep-release-105.0 && \
    perl INSTALL.pl -a a -s homo_sapiens -y GRCh38 --NO_HTSLIB --NO_UPDATE
ENV PATH="/ensembl-vep-release-105.0:${PATH}"

# Install MultiQC
RUN pip3 install multiqc

# Set working directory
WORKDIR /data

CMD ["bash"]
