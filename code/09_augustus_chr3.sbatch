#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH -J augustus_chr3
#SBATCH --output=%x.%j.out

module purge
module load AUGUSTUS/3.5.0-gfbf-2024a

export AUGUSTUS_CONFIG_PATH=$HOME/augustus_config

GENOME="analyses/02_assembly/chr3/polishing/polished_chr3.fasta"
OUTDIR="analyses/03_annotation/augustus_chr3"
OUTFILE="$OUTDIR/augustus_chr3.gff"

mkdir -p "$OUTDIR"

augustus \
  --species=arabidopsis \
  --gff3=on \
  --softmasking=1 \
  "$GENOME" \
  > "$OUTFILE"

echo "=== DONE ==="
