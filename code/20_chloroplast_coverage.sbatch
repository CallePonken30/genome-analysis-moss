#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 06:00:00
#SBATCH -J cp_cov
#SBATCH --output=logs/slurm/cp_cov.%j.out

set -euo pipefail

module load Bowtie2
module load SAMtools

echo "=== START chloroplast coverage ==="
date

PROJECT="$HOME/genome-analysis-moss"
OUTDIR="/proj/uppmax2026-1-61/nobackup/$USER/chloroplast_coverage"

REF="$PROJECT/analyses/09_organelle/chloroplast_retry/embplant_pt.K105.scaffolds.graph1.1.path_sequence.fasta"
R1="$PROJECT/data/raw_data/whole_genome/illumina/CRR809859_f1.fq.gz"
R2="$PROJECT/data/raw_data/whole_genome/illumina/CRR809859_r2.fq.gz"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

rm -f cp.bam cp.bam.bai cp.depth cp_index*

echo "=== Checking input files ==="
ls -lh "$REF" "$R1" "$R2"

echo "=== Build Bowtie2 index ==="
bowtie2-build "$REF" cp_index

echo "=== Mapping and sorting directly to BAM ==="
bowtie2 -p 4 -x cp_index \
  -1 "$R1" \
  -2 "$R2" \
| samtools sort -@ 4 -m 2G -o cp.bam

echo "=== Index BAM ==="
samtools index cp.bam

echo "=== Calculate coverage depth ==="
samtools depth cp.bam > cp.depth

echo "=== Average coverage ==="
awk '{sum+=$3} END {print "Average coverage:", sum/NR}' cp.depth

echo "=== Output files ==="
ls -lh cp.bam cp.bam.bai cp.depth

echo "=== DONE ==="
date
