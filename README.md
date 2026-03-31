# Genome Analysis Project – Paper 2

## Project
This repository is for the Genome Analysis lab project based on:

**Zhou et al. 2023**  
*Chromosome-level genome assembly of Niphotrichum japonicum provides new insights into heat stress responses in mosses.*

## Main aim
The aim is to reproduce and understand parts of the published bioinformatics workflow for the moss *Niphotrichum japonicum*, with focus on genome assembly, annotation, and differential gene expression under heat stress.

## Planned workflow
1. Quality control of sequencing reads
2. Read trimming
3. Nanopore genome assembly
4. Assembly polishing with Illumina reads
5. Assembly evaluation
6. Structural and functional annotation
7. RNA-seq mapping and differential expression analysis
8. Biological interpretation

## Repository structure

- `code/` – scripts and commands  
- `analyses/` – outputs from each analysis step  
- `data/` – metadata and symbolic links to raw data  
- `docs/` – project planning and notes  
- `figures/` – plots and images for interpretation and presentation

## Notes
Large raw sequencing data files are stored on UPPMAX and accessed via symbolic links.  
Only scripts, metadata, and summary results are stored in this repository.
