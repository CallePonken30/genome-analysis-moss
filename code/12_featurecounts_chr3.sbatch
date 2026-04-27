#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 2
#SBATCH -t 02:00:00
#SBATCH -J featurecounts_chr3
#SBATCH --output=logs/slurm/%x.%j.out

module purge
module load Subread/2.1.1-GCC-13.3.0

# Use BRAKER annotation instead of AUGUSTUS
ANNOT="analyses/03_annotation/braker3_chr3/braker_output/braker.gtf"
OUT="analyses/06_counts/featurecounts_braker_chr3.txt"

featureCounts \
  -T 2 \
  -p \
  -F GTF \
  -t exon \
  -g gene_id \
  -a "$ANNOT" \
  -o "$OUT" \
  analyses/05_expression/mapping/control/Control_1.sorted.bam \
  analyses/05_expression/mapping/control/Control_2.sorted.bam \
  analyses/05_expression/mapping/control/Control_3.sorted.bam \
  analyses/05_expression/mapping/heat_treated/Heat_1.sorted.bam \
  analyses/05_expression/mapping/heat_treated/Heat_2.sorted.bam \
  analyses/05_expression/mapping/heat_treated/Heat_3.sorted.bam

echo "=== DONE ==="
