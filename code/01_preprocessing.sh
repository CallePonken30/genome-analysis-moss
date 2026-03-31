#!/bin/bash

# ============================================
# 01_preprocessing.sh
# Quality control of raw Illumina reads
# ============================================

set -e

echo "Loading FastQC..."
module load FastQC

echo "Running FastQC on raw reads..."
fastqc -t 2 data/raw_data/illumina/*.fq.gz \
-o analyses/01_preprocessing/fastqc_raw/

echo "Raw read FastQC completed."
