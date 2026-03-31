#!/bin/bash

# ============================================
# 01_preprocessing.sh
# Quality control of raw reads with FastQC
# ============================================

set -e

echo "Loading FastQC..."
module load FastQC

# -------- WHOLE GENOME ILLUMINA --------
echo "Running FastQC on whole-genome Illumina reads..."
fastqc -t 2 data/raw_data/whole_genome/illumina/*.fq.gz \
-o analyses/01_preprocessing/fastqc_raw/whole_genome/illumina/

# -------- CHR3 ILLUMINA --------
echo "Running FastQC on chr3 Illumina reads..."
fastqc -t 2 data/raw_data/chr3/illumina/*.fastq.gz \
-o analyses/01_preprocessing/fastqc_raw/chr3/illumina/

# -------- CHR3 Hi-C --------
echo "Running FastQC on chr3 Hi-C reads..."
fastqc -t 2 data/raw_data/chr3/hic/*.fastq.gz \
-o analyses/01_preprocessing/fastqc_raw/chr3/hic/

# -------- CHR3 NANOPORE --------
echo "Running FastQC on chr3 Nanopore reads..."
fastqc -t 2 data/raw_data/chr3/nanopore/*.fq.gz \
-o analyses/01_preprocessing/fastqc_raw/chr3/nanopore/

echo "Raw read FastQC completed."
