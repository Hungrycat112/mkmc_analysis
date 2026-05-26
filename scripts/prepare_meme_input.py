import pandas as pd

"""
Prepare MEME input FASTA files from MKMC k-mer correlation output. Note: this script was also used to prepare weblogo fasta format

This script:
1. Loads MKMC top k-mer correlation results
2. Selects the strongest positively and negatively age-correlated k-mers
3. Saves selected k-mers as plain text files
4. Converts selected k-mers into FASTA format for MEME motif analysis

-> afterwards, manually padded sequences (MEME requires certain length) via this command-line: 
awk '{print ">kmer" NR "\nN" $0 "N"}' blood_top100_trimmed_neg.txt > blood_trimmed_neg_padded.fa 


"""

# Input and output settings
input_file = "1k_size_7_raw_cs_spearman_top"
top_n = 1000 # can change how many top k-mers you want to select

positive_txt = "blood_top1000_trimmed_pos.txt"
negative_txt = "blood_top1000_trimmed_neg.txt"

positive_fasta = "blood_trimmed_pos.fa"
negative_fasta = "blood_trimmed_neg.fa"


# Load MKMC correlation output
# The MKMC output is whitespace-delimited.
# Expected columns include "k-mer" and "spearman".
df = pd.read_csv(input_file, sep=r"\s+")

required_cols = {"k-mer", "spearman"}
missing_cols = required_cols - set(df.columns)

if missing_cols:
    raise ValueError(f"Missing required columns: {missing_cols}")

print(f"Loaded {len(df):,} k-mers from {input_file}")

# Select top correlated k-mers
# Positive k-mers increase in frequency with age.
top_pos = df.sort_values("spearman", ascending=False).head(top_n)

# Negative k-mers decrease in frequency with age.
top_neg = df.sort_values("spearman", ascending=True).head(top_n)

# Save selected k-mers as text
top_pos["k-mer"].to_csv(positive_txt, index=False, header=False)
top_neg["k-mer"].to_csv(negative_txt, index=False, header=False)

print(f"Saved positive k-mers: {positive_txt}")
print(f"Saved negative k-mers: {negative_txt}")

# Convert k-mer text files to FASTA
def kmers_to_fasta(input_file, output_file, prefix):
    """
    Convert one-k-mer-per-line text file into FASTA format.

    MEME requires FASTA input, where each sequence has:
    - a header line beginning with ">"
    - a sequence line containing the k-mer
    """
    with open(input_file, "r") as infile:
        kmers = infile.read().splitlines()

    with open(output_file, "w") as fasta:
        for i, kmer in enumerate(kmers, start=1):
            fasta.write(f">{prefix}_{i}\n{kmer}\n")

    print(f"Saved FASTA: {output_file}")


kmers_to_fasta(positive_txt, positive_fasta, prefix="pos_kmer")
kmers_to_fasta(negative_txt, negative_fasta, prefix="neg_kmer")
