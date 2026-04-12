#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -J polish_chr3
#SBATCH --mail-type=ALL
#SBATCH --output=analyses/02_assembly/chr3/polishing/%x.%j.out

set -euo pipefail

module purge
module load BWA/0.7.19-GCCcore-13.3.0
module load SAMtools/1.22.1-GCC-13.3.0
module load Pilon/1.24-Java-17

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
ASSEMBLY="${PROJECT_DIR}/analyses/02_assembly/chr3/assembly.fasta"
POLISH_DIR="${PROJECT_DIR}/analyses/02_assembly/chr3/polishing"
R1="${PROJECT_DIR}/data/trimmed_data/chr3/illumina/R1_paired.fastq.gz"
R2="${PROJECT_DIR}/data/trimmed_data/chr3/illumina/R2_paired.fastq.gz"
BAM="${POLISH_DIR}/chr3_aln.sorted.bam"

mkdir -p "${POLISH_DIR}"

echo "=== START POLISHING ==="
echo "Project dir: ${PROJECT_DIR}"
echo "Assembly: ${ASSEMBLY}"
echo "R1: ${R1}"
echo "R2: ${R2}"

rm -f "${POLISH_DIR}/chr3_aln.sorted.bam"*
rm -f "${POLISH_DIR}/polished_chr3"*
rm -f "${POLISH_DIR}"/*.vcf
rm -f "${POLISH_DIR}"/*.changes

test -f "${ASSEMBLY}"
test -f "${R1}"
test -f "${R2}"

bwa index "${ASSEMBLY}"

bwa mem -t 4 "${ASSEMBLY}" "${R1}" "${R2}" | \
  samtools sort -@ 4 -o "${BAM}"

test -f "${BAM}"

samtools index "${BAM}"

java -Xmx8G -jar "$EBROOTPILON/pilon.jar" \
  --genome "${ASSEMBLY}" \
  --frags "${BAM}" \
  --output polished_chr3 \
  --outdir "${POLISH_DIR}" \
  --threads 4

echo "=== DONE ==="
