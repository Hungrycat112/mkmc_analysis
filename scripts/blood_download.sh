#!/bin/bash
# This script downloads blood samples from SRA -> converts -> compresses -> validates FASTQ files
# Usage: qsub blood_download.sh
# Input: tab-separated file with SRA accession ID (one per line)
# Output: blood_samples/${sample}_R1.fastq.gz

#$ -cwd
#$ -V  
#$ -N blood_dl # job name
#$ -l h_rt=24:00:00,h_data=64G # 24hr runtime, 64 GB memory per slot
#$ -pe shared 2
#$ -t 1-5 # task array (one task per sample)
#$ -tc 5 # max is 5 concurrent tasks
#$ -o logs/blood_dl_$TASK_ID.out
#$ -e logs/blood_dl_$TASK_ID.err

export PATH=/u/scratch/s/sonya/sratoolkit.3.2.0-centos_linux64/bin:$PATH # path to mkmc
mkdir -p blood_samples logs sra_tmp # temp directories

INPUT="100_blood_input.txt" # blood input file

# get sample ID for task
sample=$(awk -v line=$SGE_TASK_ID 'NR==line {print $1}' "$INPUT")

out="blood_samples/${sample}_R1.fastq.gz" # final output path

tmpdir="fastq_tmp/${sample}_$$" # fasterq-dump output goes here
fqtmp="fasterq_tmp/${sample}_$$" # fasterq-dump internal temp files go here
mkdir -p "$tmpdir" "$fqtmp"

# exits if no sample found in task ID
if [[ -z "$sample" ]]; then
    echo "No sample for task $SGE_TASK_ID"
    exit 0
fi

# skip if output already exists + passes gzip validation
if [[ -f "$out" ]] && gzip -t "$out" 2>/dev/null; then
    echo "Skipping $sample — already done"
    exit 0
fi

echo "Processing $sample"
date

# clean up temp dirs on exit, even if job is killed
trap 'rm -rf "$tmpdir" "$fqtmp"' EXIT

if [[ ! -f "sra_tmp/$sample/$sample.sra" ]]; then
    echo "Prefetching $sample"
    rm -f "sra_tmp/$sample/$sample.sra.lock" # removes stale lock/tmp files from previous failed/killed downloads
    rm -f "sra_tmp/$sample/$sample.sra.tmp"
    prefetch "$sample" --output-directory sra_tmp --max-size 100G || exit 1
else
    echo "Using existing SRA: sra_tmp/$sample/$sample.sra"
fi

# verify SRA file exists, before proceeding
if [[ ! -f "sra_tmp/$sample/$sample.sra" ]]; then
    echo "ERROR: SRA file not found for $sample"
    exit 1
fi

# convert SRA to FASTQ
echo "Running fasterq-dump for $sample"
date
# try fasterq-dump first because faster
fasterq-dump "sra_tmp/$sample/$sample.sra" \
    --outdir "$tmpdir" \
    --temp   "$fqtmp" \
    --threads 1 \
    --mem 10000 \
    --split-files

# fall back to fastq-dump if fasterq-dump fails (lower memory usage)
if [[ $? -ne 0 ]]; then
    echo "fasterq-dump failed, trying fastq-dump fallback for $sample"
    rm -rf "$tmpdir"/*
    fastq-dump "sra_tmp/$sample/$sample.sra" \
        --outdir "$tmpdir" \
        --split-files || exit 1
fi

# compress FASTQ to gzip
echo "Compressing $sample"
date
# handles both paired-end (_1.fastq) and single-end (.fastq) output
if [[ -f "$tmpdir/${sample}_1.fastq" ]]; then
    gzip -c "$tmpdir/${sample}_1.fastq" > "$out"
elif [[ -f "$tmpdir/${sample}.fastq" ]]; then
    gzip -c "$tmpdir/${sample}.fastq" > "$out"
else
    echo "ERROR: no FASTQ found for $sample"
    ls -lh "$tmpdir"
    exit 1
fi

# validate output
gzip -t "$out" || {
    echo "ERROR: gzip validation failed for $out"
    rm -f "$out"
    exit 1
}

echo "Finished $sample"
date
# Temp dirs are automatically cleaned up by trap EXIT
