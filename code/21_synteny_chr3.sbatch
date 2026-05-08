#!/bin/bash -l
#SBATCH -A uppmax2026-1-61
#SBATCH -p pelle
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -J synteny_chr3
#SBATCH --output=logs/slurm/synteny_chr3.%j.out

set -euo pipefail

module load minimap2
module load R/4.4.2-gfbf-2024a

echo "=== START synteny analysis ==="
date

PROJECT="$HOME/genome-analysis-moss"

REF="$PROJECT/data/reference/GCA_014871385.1_CpurpureusR40_1_0_genomic.fna"

QUERY="$PROJECT/analyses/02_assembly/chr3/polishing/polished_chr3.fasta"

OUTDIR="$PROJECT/analyses/10_synteny"

mkdir -p "$OUTDIR"

echo "=== Running minimap2 alignment ==="

minimap2 -x asm10 \
"$REF" \
"$QUERY" \
> "$OUTDIR/chr3_vs_cpurpureus.paf"

echo "=== Creating synteny plot ==="

cat > "$OUTDIR/plot_synteny.R" << 'RSCRIPT'
paf <- read.table(
  "chr3_vs_cpurpureus.paf",
  sep="\t",
  header=FALSE
)

pdf(
  "chr3_vs_cpurpureus.pdf",
  width=8,
  height=8
)

plot(
  paf$V8,
  paf$V3,
  pch=16,
  cex=0.3,
  xlab="Ceratodon purpureus genome",
  ylab="N. japonicum chr3 assembly",
  main="Synteny analysis"
)

dev.off()
RSCRIPT

cd "$OUTDIR"

Rscript plot_synteny.R

echo "=== Output files ==="
ls -lh

echo "=== DONE ==="
date
