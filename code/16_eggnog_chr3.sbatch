#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 2
#SBATCH -t 24:00:00
#SBATCH -J eggnog_braker
#SBATCH --output=logs/slurm/%x.%j.out

module purge
module load eggnog-mapper/2.1.13-gfbf-2024a

# Input from BRAKER3
INPUT_PROTEINS="/home/calle30/genome-analysis-moss/analyses/03_annotation/braker3_chr3/braker_output/braker.aa"

# Output folder
OUT_DIR="/home/calle30/genome-analysis-moss/analyses/03_annotation/eggnog_braker_chr3"

# Shared UPPMAX eggNOG database
DATA_DIR="/sw/data/eggNOG/5.0.0/rackham/"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

emapper.py \
  -i "$INPUT_PROTEINS" \
  --output chr3_braker_functional \
  --cpu 2 \
  --itype proteins \
  --data_dir "$DATA_DIR"

echo "=== DONE ==="
