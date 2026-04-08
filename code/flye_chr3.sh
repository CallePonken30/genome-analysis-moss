#!/bin/bash
#SBATCH -A uppmax2026-1-61
#SBATCH -J flye_chr3
#SBATCH -t 12:00:00
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -o slurm-flye_chr3.out

module load Flye/2.9.6-GCC-13.3.0

flye \
  --nano-raw data/raw_data/chr3/nanopore/chr3_clean_nanopore.fq.gz \
  --out-dir analyses/02_assembly/chr3 \
  --threads 4 \
  --genome-size 500m \
  --asm-coverage 40
