#!/bin/bash
#$ -cwd
#$ -V
#$ -N trim_1k
#$ -l h_rt=24:00:00,h_data=6G
#$ -pe shared 4
#$ -o logs/trim_1k.out
#$ -e logs/trim_1k.err

# trim_1k.sh
#
# Two-step trimming of dog FASTQ files:
#   1. fastp quality filtering (Phred >= 20, min length 7bp)
#   2. Hard truncation to first 7 bases via awk
#
# Input:  1k_input.txt (col 2 = filenames), 1k_dogs/*.fastq.gz
# Output: trimmed_fastq/*_7bp.fastq.gz, logs/
#
# Safe to rerun - skips already completed samples
#
# Usage: qsub trim_1k.sh
# Dependencies: fastp, awk, gzip


mkdir -p trimmed_fastq
mkdir -p temp_trimmed
mkdir -p logs

COUNT=0
BATCH_SIZE=10

while read fq; do

    input="1k_dogs/$fq"

    base=$(basename "$fq" .fastq.gz)

    temp_out="temp_trimmed/${base}_fastp.fastq.gz"
    final_out="trimmed_fastq/${base}_7bp.fastq.gz"

    # Skip if already processed successfully
    if [[ -f "$final_out" ]] && gzip -t "$final_out" 2>/dev/null; then
        echo "Skipping $fq (already done)"
        rm -f "$input"
        COUNT=$((COUNT+1))
        continue
    fi

    # Check input exists
    if [[ ! -f "$input" ]]; then
        echo "Missing: $input"
        continue
    fi

    echo "Processing $fq"

    # Step 1: fastp quality trimming
    fastp \
        -i "$input" \
        -o "$temp_out" \
        --disable_adapter_trimming \
        --thread 4 \
        --qualified_quality_phred 20 \
        --length_required 7 \
        --json logs/${base}.json \
        --html logs/${base}.html

    # Check fastp succeeded
    if [[ ! -s "$temp_out" ]]; then
        echo "ERROR: fastp failed on $fq"
        continue
    fi

    # Step 2: keep only first 7 bases
    zcat "$temp_out" | awk '
        NR%4==1 {print}
        NR%4==2 {print substr($0,1,7)}
        NR%4==3 {print}
        NR%4==0 {print substr($0,1,7)}
    ' | gzip > "$final_out"

    # Validate final output before deleting originals, delete originals to maintain storage
    if [[ -s "$final_out" ]] && gzip -t "$final_out" 2>/dev/null; then

        rm -f "$temp_out"
        rm -f "$input"

        echo "Finished $fq"

    else

        echo "ERROR: invalid output for $fq"
        rm -f "$final_out"

    fi

    COUNT=$((COUNT+1))

    if (( COUNT % BATCH_SIZE == 0 )); then
        echo "Completed $COUNT samples"
    fi

done < <(awk '{print $2}' 1k_input.txt)

echo "All done. Total processed: $COUNT"
