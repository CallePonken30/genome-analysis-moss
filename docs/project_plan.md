# Project plan – Genome Analysis (Paper 2)

## Project overview
This project aims to reproduce and understand key parts of the bioinformatics analyses presented in the paper:

Zhou et al. (2023) – Chromosome-level genome assembly of *Niphotrichum japonicum* provides new insights into heat stress responses in mosses.

The focus of the project is to reconstruct the genome assembly workflow, evaluate genome quality, perform genome annotation, and analyse gene expression under heat stress conditions.

---

## Goals
The main goals of the project are:

- Learn how genome assembly workflows are performed in practice
- Evaluate the quality and completeness of a genome assembly
- Perform structural and functional genome annotation
- Analyse RNA-seq data to identify genes responding to heat stress
- Interpret biological results in an evolutionary and ecological context
- Develop skills in reproducible bioinformatics using GitHub and UPPMAX

---

## Type of samples
The study organism is the moss *Niphotrichum japonicum*.

Samples include:

- DNA sequencing reads used for genome assembly
- RNA-seq reads collected from moss exposed to heat stress and control conditions

---

## Type of data
The project will use several sequencing data types:

- Nanopore long-read sequencing data for primary genome assembly
- Illumina short-read sequencing data for assembly polishing
- RNA-seq data for differential gene expression analysis
- Possibly Hi-C data for chromosome-level scaffolding

---

## Planned analyses and software
The planned workflow includes:

1. Quality control of sequencing reads  
   Software: FastQC  

2. Read trimming and filtering  
   Software: Trimmomatic or similar  

3. Genome assembly of long reads  
   Software: Flye or Canu  

4. Assembly polishing with short reads  
   Software: Pilon  

5. Assembly quality evaluation  
   Software: BUSCO and QUAST  

6. Structural genome annotation  
   Software: BRAKER / AUGUSTUS  

7. Functional annotation  
   Software: InterProScan or similar  

8. RNA-seq mapping  
   Software: HISAT2 or STAR  

9. Differential expression analysis  
   Software: DESeq2  

---

## Time frame and checkpoints
The project will follow the course timeline:

- Early April: repository setup and project planning
- Mid April: genome assembly and polishing
- Late April: assembly evaluation and annotation
- Early May: RNA-seq analysis and biological interpretation
- Mid May: preparation of presentation and final documentation

Important checkpoints include:

- Successful genome assembly
- Assembly quality metrics obtained
- Annotation files generated
- Differential expression results produced

---

## Plan for long running analyses
Genome assembly and RNA-seq mapping may require significant computational time.

To manage this:

- Jobs will be submitted to the UPPMAX cluster using SLURM
- Scripts will be tested interactively before batch submission
- Log files will be stored for troubleshooting
- Intermediate outputs will be organised in analysis-specific folders

---

## Data management plan
Raw sequencing data will remain stored in the shared project directory on UPPMAX.

Only symbolic links will be created in the working directory to avoid unnecessary data duplication.

Large intermediate files will be stored in the project work directory rather than the GitHub repository.

The GitHub repository will contain:

- scripts
- documentation
- metadata
- summary outputs
- figures

---

## Expected results
It is expected that:

- a draft genome assembly of reasonable completeness will be produced
- a set of predicted genes will be identified
- a subset of genes will show differential expression under heat stress
- results will be broadly consistent with the conclusions of the original paper, although minor differences may occur due to methodological choices

---

## Potential challenges
Possible challenges include:

- large computational resource requirements
- software installation and compatibility issues
- long runtimes of genome assembly
- interpretation of complex biological results

These challenges will be addressed through documentation, collaboration with peers, and consultation with teaching assistants during lab sessions.
