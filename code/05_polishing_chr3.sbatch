#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -J polish_chr3
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out

module purge
module spider bwa
module spider samtools
module spider pilon

# Index assembly
bwa index analyses/02_assembly/chr3/assembly.fasta

# Map trimmed Illumina reads to assembly and sort BAM
bwa mem -t 4 \
  analyses/02_assembly/chr3/assembly.fasta \
  data/trimmed_data/chr3/illumina/R1_paired.fastq.gz \
  data/trimmed_data/chr3/illumina/R2_paired.fastq.gz | \
  samtools sort -@ 4 -o analyses/02_assembly/chr3/polishing/chr3_aln.sorted.bam

# Index BAM
samtools index analyses/02_assembly/chr3/polishing/chr3_aln.sorted.bam

# Run Pilon
pilon \
  --genome analyses/02_assembly/chr3/assembly.fasta \
  --frags analyses/02_assembly/chr3/polishing/chr3_aln.sorted.bam \
  --output polished_chr3 \
  --outdir analyses/02_assembly/chr3/polishing \
  --threads 4
