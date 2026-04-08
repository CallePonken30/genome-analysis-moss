#!/bin/bash

module load Flye

mkdir -p analyses/02_assembly/chr3

flye \
  --nano-raw data/raw_data/chr3/nanopore/chr3_clean_nanopore.fq.gz \
  --out-dir analyses/02_assembly/chr3 \
  --threads 2
