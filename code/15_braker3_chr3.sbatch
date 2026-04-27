#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 15:00:00
#SBATCH --mem=64G
#SBATCH -J braker3_chr3
#SBATCH --output=logs/slurm/%x.%j.out

set -e
# -------------------------
# File variables
# -------------------------

OUT_DIR="/home/calle30/genome-analysis-moss/analyses/03_annotation/braker3_chr3"

GENOME="/home/calle30/genome-analysis-moss/analyses/02_assembly/chr3/polishing/polished_chr3.fasta"

BAM_FILE="/home/calle30/genome-analysis-moss/analyses/05_expression/mapping/merged/rna_all.sorted.bam"

BRAKER_SIF="/proj/uppmax2026-1-61/Genome_Analysis/2_Zhou_2023/braker3.sif"

AUGUSTUS_CONFIG="/home/calle30/augustus_config"

# -------------------------
# Setup
# -------------------------

mkdir -p "$OUT_DIR"
mkdir -p "$OUT_DIR/tmp"
mkdir -p logs/slurm

export TMPDIR="$OUT_DIR/tmp"

# -------------------------
# Run BRAKER3
# -------------------------

singularity exec \
  -B /home/calle30:/home/calle30 \
  -B "$AUGUSTUS_CONFIG":/opt/Augustus/config \
  "$BRAKER_SIF" \
  braker.pl \
    --genome="$GENOME" \
    --bam="$BAM_FILE" \
    --species=calle_chr3_001 \
    --threads=4 \
    --min_contig=5000 \
    --skipOptimize \
    --workingdir="$OUT_DIR/braker_output"

echo "=== DONE ==="
