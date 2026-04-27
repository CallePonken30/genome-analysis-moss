#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 06:00:00
#SBATCH --mem=32G
#SBATCH -J getorg_cp
#SBATCH --output=logs/slurm/%x.%j.out

module purge
module load GetOrganelle/1.7.7.1-foss-2024a

R1="data/raw_data/whole_genome/illumina/CRR809859_f1.fq.gz"
R2="data/raw_data/whole_genome/illumina/CRR809859_r2.fq.gz"

OUTDIR="analyses/09_organelle/chloroplast_retry"

get_organelle_from_reads.py \
  -1 "$R1" \
  -2 "$R2" \
  -o "$OUTDIR" \
  -F embplant_pt \
  -R 25 \
  -k 21,45,65,85,105 \
  --reduce-reads-for-coverage 200 \
  -t 4

echo "=== DONE chloroplast retry ==="
