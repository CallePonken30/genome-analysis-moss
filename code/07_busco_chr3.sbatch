#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 2
#SBATCH -t 04:00:00
#SBATCH -J busco_chr3
#SBATCH --output=%x.%j.out

module load BUSCO/5.8.2-gfbf-2024a

ASSEMBLY="analyses/02_assembly/chr3/polishing/polished_chr3.fasta"
BUSCO_DIR="analyses/04_evaluation/busco_chr3"
LINEAGE="embryophyta_odb10"

mkdir -p "${BUSCO_DIR}"

test -f "${ASSEMBLY}"

busco \
  -i "${ASSEMBLY}" \
  -o busco_chr3 \
  -m genome \
  -l "${LINEAGE}" \
  -c 2 \
  --out_path "${BUSCO_DIR}"

test -d "${BUSCO_DIR}/busco_chr3"

echo "=== DONE ==="
