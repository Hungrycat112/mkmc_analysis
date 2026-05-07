# mkmc_analysis_pipeline

## Overview
This repository has computaional workflows for identifying age-associated k-mer patterns in canine sequencing data, buccal and blood, using mkmc with machine learning analysis. 

## Pipeline overview: 
FASTQ reads -> preprocess / trim -> mkmc -> spearman -> top k-mer selected -> machine learning models -> evaluation + visualization

## Methods: 
typical mkmc command: ./mkmc -m10 -k 7 -t 12 --cs 1000000000 --cor spearman -n freq -o matrix -f fq -p phenotype.txt input.txt output_prefix temp/

## Input Format
### input.txt

Tab-delimited file containing:
1. Sample identifier
2. Path to FASTQ file

Example:

```text
31231210401824    1k_dogs/DW31231210401824_R1.fastq.gz
31231210401827    1k_dogs/DW31231210401827_R1.fastq.gz
31240210300924    1k_dogs/DT31240210300924_R1.fastq.gz
```

### phenotype.txt

One chronological age per line corresponding to the sample order in `input.txt`.

Example:

```text
2
7
10
```

##  Machine learning: 
Ridge Regression model + LOOCV

## Motif Analysis: 
- weblogo
- logomaker
- MEME suite


