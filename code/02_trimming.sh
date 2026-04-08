#!/bin/bash
#SBATCH -A uppmax2026-1-61
#SBATCH --qos=short
#SBATCH -t 01:00:00
#SBATCH -c 2
#SBATCH -J trimming_all

module load Trimmomatic

echo "=== START TRIMMING ==="

# =========================
# CHR3 ILLUMINA
# =========================
echo "Trimming chr3 illumina..."

trimmomatic PE -threads 2 \
data/raw_data/chr3/illumina/chr3_illumina_R1.fastq.gz \
data/raw_data/chr3/illumina/chr3_illumina_R2.fastq.gz \
data/trimmed_data/chr3/illumina/R1_paired.fastq.gz \
data/trimmed_data/chr3/illumina/R1_unpaired.fastq.gz \
data/trimmed_data/chr3/illumina/R2_paired.fastq.gz \
data/trimmed_data/chr3/illumina/R2_unpaired.fastq.gz \
SLIDINGWINDOW:4:20 MINLEN:50

# =========================
# WHOLE GENOME ILLUMINA
# =========================
echo "Trimming whole genome illumina..."

trimmomatic PE -threads 2 \
data/raw_data/whole_genome/illumina/CRR809859_f1.fq.gz \
data/raw_data/whole_genome/illumina/CRR809859_r2.fq.gz \
data/trimmed_data/whole_genome/illumina/R1_paired.fastq.gz \
data/trimmed_data/whole_genome/illumina/R1_unpaired.fastq.gz \
data/trimmed_data/whole_genome/illumina/R2_paired.fastq.gz \
data/trimmed_data/whole_genome/illumina/R2_unpaired.fastq.gz \
SLIDINGWINDOW:4:20 MINLEN:50

echo "=== DONE ==="
