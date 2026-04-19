#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 12:00:00
#SBATCH -J braker_chr3
#SBATCH --output=%x.%j.out

module purge
module load AUGUSTUS/3.5.0-gfbf-2024a
module load GeneMark-ET/4.72-GCCcore-13.3.0
module load BRAKER/3.0.8-foss-2024a

export AUGUSTUS_CONFIG_PATH=$HOME/augustus_config

GENOME="analyses/02_assembly/chr3/polishing/polished_chr3.fasta"
BAM="analyses/05_expression/mapping/merged/rna_all.sorted.bam"
OUTDIR="analyses/03_annotation/braker_chr3"

# Create output directory
mkdir -p "$OUTDIR"

braker.pl \
  --genome="$GENOME" \
  --bam="$BAM" \
  --species=chr3_moss \
  --softmasking \
  --threads=4 \
  --workingdir="$OUTDIR"

echo "=== DONE ==="
