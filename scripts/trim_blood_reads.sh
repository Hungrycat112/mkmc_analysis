#!/bin/bash
#$ -cwd
#$ -V
#$ -N trim_blood
#$ -l h_rt=24:00:00,h_data=6G
#$ -pe shared 16
#$ -o logs/trim_blood.out
#$ -e logs/trim_blood.err
# trim_blood.sh
#
# Two-step trimming of blood FASTQ files:
#   1. fastp quality filtering (Phred >= 20, min length 7bp)
#   2. Hard truncation to first 7 bases via awk
#
# Input:  blood_samples/*.fastq.gz
# Output: trimmed_blood/*_7bp.fastq.gz, logs/
#
# Safe to rerun - skips already completed samples
# Processes PARALLEL_JOBS samples simultaneously
#
# Usage: qsub trim_blood.sh
# Dependencies: fastp, awk, gzip

mkdir -p trimmed_blood # final trimmed finals go here
mkdir -p temp_trimmed_blood
mkdir -p logs

PARALLEL_JOBS=4       # number of samples to process at once
THREADS_PER_JOB=4     # fastp threads per sample (PARALLEL_JOBS * THREADS_PER_JOB should = pe shared above)
COUNT=0
BATCH_SIZE=10

process_sample() {
    input="$1"
    fq=$(basename "$input")
    base=$(basename "$fq" .fastq.gz)
    temp_out="temp_trimmed_blood/${base}_fastp.fastq.gz"
    final_out="trimmed_blood/${base}_7bp.fastq.gz"

    # Skip if already processed successfully
    if [[ -f "$final_out" ]] && gzip -t "$final_out" 2>/dev/null; then
        echo "Skipping $fq (already done)"
        return
    fi

    # Check input exists
    if [[ ! -f "$input" ]]; then
        echo "Missing: $input"
        return
    fi

    echo "Processing $fq"

    # Step 1: fastp quality trimming
    fastp \
        -i "$input" \
        -o "$temp_out" \
        --disable_adapter_trimming \
        --thread "$THREADS_PER_JOB" \
        --qualified_quality_phred 20 \
        --length_required 7 \
        --json logs/${base}.json \
        --html logs/${base}.html

    # Check fastp succeeded
    if [[ ! -s "$temp_out" ]]; then
        echo "ERROR: fastp failed on $fq"
        return
    fi

    # Step 2: keep only first 7 bases
    zcat "$temp_out" | awk '
        NR%4==1 {print}
        NR%4==2 {print substr($0,1,7)}
        NR%4==3 {print}
        NR%4==0 {print substr($0,1,7)}
    ' | gzip > "$final_out"

    # Validate final output
    if [[ -s "$final_out" ]] && gzip -t "$final_out" 2>/dev/null; then
        rm -f "$temp_out"
        echo "Finished $fq"
    else
        echo "ERROR: invalid output for $fq"
        rm -f "$final_out"
    fi
}

export -f process_sample

# Run in parallel using background jobs
JOBS_RUNNING=0

for input in blood_samples/*_R1.fastq.gz; do
    process_sample "$input" &
    JOBS_RUNNING=$((JOBS_RUNNING + 1))
    COUNT=$((COUNT + 1))

    # Wait when we hit the parallel job limit
    if (( JOBS_RUNNING >= PARALLEL_JOBS )); then
        wait
        JOBS_RUNNING=0
    fi

    if (( COUNT % BATCH_SIZE == 0 )); then
        echo "Launched $COUNT samples"
    fi
done

# Wait for any remaining background jobs to finish
wait

echo "All done. Total processed: $COUNT"
