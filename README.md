# mkmc_analysis_pipeline

## Overview
This repository has computaional workflows for identifying age-associated k-mer patterns in canine sequencing data, buccal and blood, using mkmc with machine learning analysis. 

## Pipeline overview: 
FASTQ reads -> preprocess / trim -> mkmc -> spearman -> top k-mer selected -> machine learning models -> evaluation + visualization

## Methods: 
typical mkmc command: ./mkmc -m10 -k 7 -t 12 --cs 1000000000 --cor spearman -n freq -o matrix -f fq -p phenotype.txt input.txt output_prefix temp/

##  Machine learning: 
Ridge Regression model + LOOCV

## Motif Analysis: 
- weblogo
- logomaker
- MEME suite


