#!/bin/bash
# this script runs mkmc 
# usage: qsub run_mkmc.sh
# input: input.txt (list of  FASTQ files), 1k_phenotype.txt (phenotype file)
# output: matrix files in current directory
#$ -cwd
#$ -V
#$ -N mkmc_k8
#$ -l h_rt=12:00:00,h_data=8G
#$ -pe shared 12
#$ -o logs/mkmc_k8.out
#$ -e logs/mkmc_k8.err

chmod +x ./mkmc

# -m10        : 10 GB memory for MKMC internal use
# -t 12       : 12 threads (matches -pe shared 12)
# -k 7        : k-mer size of 7
# --cs 1000000000 : counter size 1 billion
# -n freq     : normalize by frequency
# --cor spearman : use Spearman correlation
# -p          : phenotype file
# -o matrix   : output prefix
# -f fq       : input format is FASTQ

./mkmc -m10 -t 12 -k 7 --cs 1000000000 -n freq --cor spearman -p 1k_phenotype.txt -o matrix -f fq trimmed_input.txt 1k_size_7_trim temp/
