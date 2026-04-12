#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 2
#SBATCH -t 01:00:00
#SBATCH -J eval_chr3
#SBATCH --output=%x.%j.out

module load QUAST/5.3.0-gfbf-2024a

ASSEMBLY="analyses/02_assembly/chr3/assembly.fasta"
POLISHED="analyses/02_assembly/chr3/polishing/polished_chr3.fasta"
EVAL_DIR="analyses/04_evaluation/quast_chr3_polish"

mkdir -p "${EVAL_DIR}"

test -f "${ASSEMBLY}"
test -f "${POLISHED}"

quast.py \
  "${ASSEMBLY}" \
  "${POLISHED}" \
  -o "${EVAL_DIR}" \
  --threads 2

test -f "${EVAL_DIR}/report.txt"

echo "=== DONE ==="
