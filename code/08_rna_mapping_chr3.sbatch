#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 08:00:00
#SBATCH -J rna_map_chr3
#SBATCH --output=%x.%j.out

module load HISAT2/2.2.1-gompi-2024a
module load SAMtools/1.21-GCC-13.3.0

GENOME="analyses/02_assembly/chr3/polishing/polished_chr3.fasta"
INDEX="analyses/05_expression/hisat2_index/chr3_index"

CTRL_OUT="analyses/05_expression/mapping/control"
HEAT_OUT="analyses/05_expression/mapping/heat_treated"
MERGED_OUT="analyses/05_expression/mapping/merged"

mkdir -p "$CTRL_OUT" "$HEAT_OUT" "$MERGED_OUT"

test -f "$GENOME"

hisat2-build "$GENOME" "$INDEX"

# Control samples
hisat2 -x "$INDEX" \
  -1 data/raw_data/rna/control/Control_1_f1.fq.gz \
  -2 data/raw_data/rna/control/Control_1_r2.fq.gz \
  -p 4 | samtools sort -@ 4 -o "$CTRL_OUT/Control_1.sorted.bam"
samtools index "$CTRL_OUT/Control_1.sorted.bam"

hisat2 -x "$INDEX" \
  -1 data/raw_data/rna/control/Control_2_f1.fq.gz \
  -2 data/raw_data/rna/control/Control_2_r2.fq.gz \
  -p 4 | samtools sort -@ 4 -o "$CTRL_OUT/Control_2.sorted.bam"
samtools index "$CTRL_OUT/Control_2.sorted.bam"

hisat2 -x "$INDEX" \
  -1 data/raw_data/rna/control/Control_3_f1.fq.gz \
  -2 data/raw_data/rna/control/Control_3_r2.fq.gz \
  -p 4 | samtools sort -@ 4 -o "$CTRL_OUT/Control_3.sorted.bam"
samtools index "$CTRL_OUT/Control_3.sorted.bam"

# Heat-treated samples
hisat2 -x "$INDEX" \
  -1 data/raw_data/rna/heat_treated/Heat_treated_42_12h_1_f1.fq.gz \
  -2 data/raw_data/rna/heat_treated/Heat_treated_42_12h_1_r2.fq.gz \
  -p 4 | samtools sort -@ 4 -o "$HEAT_OUT/Heat_1.sorted.bam"
samtools index "$HEAT_OUT/Heat_1.sorted.bam"

hisat2 -x "$INDEX" \
  -1 data/raw_data/rna/heat_treated/Heat_treated_42_12h_2_f1.fq.gz \
  -2 data/raw_data/rna/heat_treated/Heat_treated_42_12h_2_r2.fq.gz \
  -p 4 | samtools sort -@ 4 -o "$HEAT_OUT/Heat_2.sorted.bam"
samtools index "$HEAT_OUT/Heat_2.sorted.bam"

hisat2 -x "$INDEX" \
  -1 data/raw_data/rna/heat_treated/Heat_treated_42_12h_3_f1.fq.gz \
  -2 data/raw_data/rna/heat_treated/Heat_treated_42_12h_3_r2.fq.gz \
  -p 4 | samtools sort -@ 4 -o "$HEAT_OUT/Heat_3.sorted.bam"
samtools index "$HEAT_OUT/Heat_3.sorted.bam"

# Merge all RNA evidence for BRAKER
samtools merge -@ 4 "$MERGED_OUT/rna_all.sorted.bam" \
  "$CTRL_OUT/Control_1.sorted.bam" \
  "$CTRL_OUT/Control_2.sorted.bam" \
  "$CTRL_OUT/Control_3.sorted.bam" \
  "$HEAT_OUT/Heat_1.sorted.bam" \
  "$HEAT_OUT/Heat_2.sorted.bam" \
  "$HEAT_OUT/Heat_3.sorted.bam"

samtools index "$MERGED_OUT/rna_all.sorted.bam"

echo "=== DONE ==="
