#!/bin/bash

module load FastQC

fastqc \
data/trimmed_data/chr3/illumina/R1_paired.fastq.gz \
data/trimmed_data/chr3/illumina/R2_paired.fastq.gz \
-o analyses/01_preprocessing/fastqc_trim/chr3/illumina
