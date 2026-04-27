#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 8
#SBATCH -t 06:00:00
#SBATCH --mem=32G
#SBATCH -J getorg_cp
#SBATCH --output=logs/slurm/%x.%j.out

module purge
module load GetOrganelle/1.7.7.1-foss-2024a

# Input reads
R1="data/raw_data/whole_genome/illumina/CRR809859_f1.fq.gz"
R2="data/raw_data/whole_genome/illumina/CRR809859_r2.fq.gz"

# Output (DO NOT pre-create this folder)
OUTDIR="analyses/09_organelle/chloroplast"

get_organelle_from_reads.py \
  -1 "$R1" \
  -2 "$R2" \
  -o "$OUTDIR" \
  -F embplant_pt \
  -R 15 \
  -k 21,45,65,85,105 \
  -t 8

echo "=== DONE chloroplast ==="
