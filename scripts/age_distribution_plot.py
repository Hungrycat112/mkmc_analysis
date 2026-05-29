import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""
Age distribution visualization for canine sequencing samples.
This script:
1. Loads phenotype files for buccal and cfDNA samples
2. Cleans age values by removing missing/non-numeric entries
3. Computes summary statistics for each tissue type
4. Generates side-by-side age distribution histograms
5. Saves publication-ready PNG and PDF figures
"""

# File paths
buccal_file = "1k_phenotype.txt"
cfdna_file = "same_blood_proj_phenotype.txt"

out_dir = "results"
os.makedirs(out_dir, exist_ok=True)

# Load phenotype files, each phenotype file has one chronological age per line
buccal = pd.read_csv(buccal_file, header=None, names=["age"])
cfdna = pd.read_csv(cfdna_file, header=None, names=["age"])

# Convert age values to numeric format
# Non-numeric values are converted to NaN and removed below
buccal["age"] = pd.to_numeric(buccal["age"], errors="coerce")
cfdna["age"] = pd.to_numeric(cfdna["age"], errors="coerce")

buccal = buccal.dropna(subset=["age"])
cfdna = cfdna.dropna(subset=["age"])

# Print sample counts + summaries
print("Buccal samples:", len(buccal))
print(buccal["age"].describe())

print("\ncfDNA samples:", len(cfdna))
print(cfdna["age"].describe())

# Summary stats helper
def summary_text(df):
    ages = df["age"]
    q1 = ages.quantile(0.25)
    q3 = ages.quantile(0.75)
    return (
        f"n = {len(ages)}\n"
        f"Median = {ages.median():.2f}\n"
        f"IQR = {q1:.2f}–{q3:.2f}\n"
        f"Range = {ages.min():.2f}–{ages.max():.2f}"
    )

buccal_summary = summary_text(buccal)
cfdna_summary = summary_text(cfdna)

# Shared bins and x-axis, use the same age bins for both tissues so distributions can be compared directly
max_age = np.ceil(max(buccal["age"].max(), cfdna["age"].max()))
bins = np.arange(0, max_age + 2, 1)

# Create side-by-side histograms
fig, axes = plt.subplots(1, 2, figsize=(8.6, 5.6), sharex=True)

# Buccal histogram
counts_buccal, _, _ = axes[0].hist(
    buccal["age"],
    bins=bins,
    edgecolor="black",
    linewidth=0.8
)

axes[0].set_title("Buccal", fontsize=13, fontweight="bold")
axes[0].set_xlabel("Age", fontsize=11)
axes[0].set_ylabel("Number of samples", fontsize=11)
axes[0].set_xlim(0, max_age + 1)
axes[0].set_ylim(0, max(counts_buccal) * 1.15)  # extra headroom
axes[0].grid(axis="y", alpha=0.2)

# cfDNA histogram
counts_cfdna, _, _ = axes[1].hist(
    cfdna["age"],
    bins=bins,
    edgecolor="black",
    linewidth=0.8
)

axes[1].set_title("Cell-free DNA (cfDNA)", fontsize=13, fontweight="bold")
axes[1].set_xlabel("Age", fontsize=11)
axes[1].set_ylabel("Number of samples", fontsize=11)
axes[1].set_xlim(0, max_age + 1)
axes[1].set_ylim(0, max(counts_cfdna) * 1.15)  # extra headroom
axes[1].grid(axis="y", alpha=0.2)

# Clean up plot borders
for ax in axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", labelsize=10)

# Adjust layout to leave room for summary boxes and caption below
plt.subplots_adjust(
    bottom=0.43,
    top=0.88,
    wspace=0.25
)

# Place summary statistic boxes under each subplot
pos0 = axes[0].get_position()
pos1 = axes[1].get_position()

fig.text(
    (pos0.x0 + pos0.x1) / 2,
    0.24,
    buccal_summary,
    ha="center",
    va="top",
    fontsize=9,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="black", linewidth=0.8)
)

fig.text(
    (pos1.x0 + pos1.x1) / 2,
    0.24,
    cfdna_summary,
    ha="center",
    va="top",
    fontsize=9,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="black", linewidth=0.8)
)

# Figure caption
fig.text(
    0.5,
    0.045,
    r"$\bf{Figure\ 2.}$ Age distributions of canine buccal and cell-free DNA (cfDNA) cohorts",
    ha="center",
    va="bottom",
    fontsize=11
)

# Save figure as PNG and PDF 
png_path = os.path.join(out_dir, "figure2_age_distribution.png")
pdf_path = os.path.join(out_dir, "figure2_age_distribution.pdf")

plt.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.15)
plt.savefig(pdf_path, bbox_inches="tight", pad_inches=0.15)
