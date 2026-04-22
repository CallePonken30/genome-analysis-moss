#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 1
#SBATCH -t 01:00:00
#SBATCH -J deseq2_chr3
#SBATCH --output=%x.%j.out

module purge
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2

Rscript code/13_deseq2_chr3.R

echo "=== DONE ==="
