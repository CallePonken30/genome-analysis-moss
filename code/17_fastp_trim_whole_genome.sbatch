#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH --mem=16G
#SBATCH -J fastp_trim
#SBATCH --output=logs/slurm/%x.%j.out

module purge
module load fastp/1.0.1-GCC-13.3.0

# Input reads
R1="data/raw_data/whole_genome/illumina/CRR809859_f1.fq.gz"
R2="data/raw_data/whole_genome/illumina/CRR809859_r2.fq.gz"

# Output directory
OUTDIR="data/trimmed_data/whole_genome/illumina"
mkdir -p $OUTDIR

OUT1="$OUTDIR/R1_paired.fastq.gz"
OUT2="$OUTDIR/R2_paired.fastq.gz"

# Report directory
REPORTDIR="analyses/01_preprocessing/fastp"
mkdir -p $REPORTDIR

fastp \
  -i $R1 \
  -I $R2 \
  -o $OUT1 \
  -O $OUT2 \
  -q 20 \
  -l 50 \
  -w 4 \
  -h $REPORTDIR/fastp_report.html \
  -j $REPORTDIR/fastp_report.json

echo "=== DONE trimming ==="
