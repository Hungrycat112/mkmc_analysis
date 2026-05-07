import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr
import os

"""
This code takes mkmc k-mer frequency matrices

1. Loads k-mer frequenc matrices from mkmc
2. Computers age-associated k-mer correlatons
3. Ridge regression with leave-one-out cross validation 
4. evaluates prediction performance
5. generates figures + summary statsitics

"""

# Global pipeline settings
USE_LOG_TRANSFORM = False
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)

# Ridge settings,  ridgeCV automatically selects the optimal alpha within each LOOCV fold
# No StandardScaler used.
RIDGE_ALPHAS = np.logspace(-4, 4, 30)

# plots
def set_pub_style():
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "legend.title_fontsize": 7,

        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "xtick.direction": "out",
        "ytick.direction": "out",

        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,

        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.transparent": False,

        "patch.linewidth": 0.5,
    })

PAL = {
    "blue": "#2166AC",
    "red": "#D6604D",
    "gray": "#888888",
    "lgray": "#CCCCCC",
    "orange": "#F4A582",
    "green": "#4DAC26",
    "bonfcol": "#B2182B",
    "p05col": "#4393C3",
}


# Scatterplot: predicted vs true
def plot_scatter(y_true, y_pred, tag="norm", out_dir="results"):
    set_pub_style()
    os.makedirs(out_dir, exist_ok=True)

    # compute regression and correlation metrics, evaulate prediction performance
    r2 = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mae = np.mean(np.abs(y_true - y_pred))
    rho, pval = spearmanr(y_true, y_pred)

    fig, ax = plt.subplots(figsize=(3.4, 3.4))

    ax.scatter(
        y_true,
        y_pred,
        s=18,
        color=PAL["blue"],
        alpha=0.65,
        linewidths=0
    )

    lims = [
        min(y_true.min(), y_pred.min()) - 0.5,
        max(y_true.max(), y_pred.max()) + 0.5
    ]

    ax.plot(
        lims,
        lims,
        "--",
        color=PAL["gray"],
        lw=0.9,
        zorder=0,
        label="Identity"
    )

    # ordinaary least square fit line used to visualize, systematic prediction bias or trend
    m, b = np.polyfit(y_true, y_pred, 1)
    ax.plot(
        np.array(lims),
        m * np.array(lims) + b,
        "-",
        color=PAL["red"],
        lw=1.2,
        zorder=1,
        label="OLS fit"
    )

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Chronological age (years)")
    ax.set_ylabel("Predicted age (years)")

    stats_text = (
        r"Spearman $\rho$" + f" = {rho:.3f}\n"
        f"$R^2$ = {r2:.3f}\n"
        f"RMSE = {rmse:.2f}\n"
        f"MAE = {mae:.2f}\n"
        f"$p$ = {pval:.2e}"
    )

    ax.text(
        0.04,
        0.97,
        stats_text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=6.2,
        bbox=dict(
            boxstyle="round,pad=0.35",
            fc="white",
            ec=PAL["lgray"],
            lw=0.6
        )
    )

    ax.legend(frameon=False, loc="lower right", handlelength=1.5)

    plt.tight_layout()

    for ext in ("pdf", "png", "svg"):
        path = os.path.join(out_dir, f"scatter_predicted_vs_true_{tag}.{ext}")
        plt.savefig(path)
        print(f"Saved: {path}")

    plt.close()

# K-mer correlation plots
def plot_kmer_correlations(
    matrix_file,
    phenotype_file,
    use_log=False,
    top_n=20,
    out_dir="results",
    tag="kmer"
):

    set_pub_style()
    os.makedirs(out_dir, exist_ok=True)

    # load mkmc matrix, rows should be samples, columns should be k-mers
    try:
        X = pd.read_csv(matrix_file, sep="\t", index_col=0)
    except Exception:
        X = pd.read_csv(matrix_file, sep=",", index_col=0)

    # some mkmc outputs are transposed depending on format, ensure samples are rows and k-mers are columns
    if X.shape[0] > X.shape[1]:
        print("Transposing k-mer matrix for correlation analysis")
        X = X.T

    y = pd.read_csv(phenotype_file, header=None).iloc[:, 0].values

    if len(y) != X.shape[0]:
        print(
            f"WARNING: sample mismatch. Matrix has {X.shape[0]} samples, "
            f"phenotype has {len(y)} ages."
        )
        min_s = min(len(y), X.shape[0])
        X = X.iloc[:min_s]
        y = y[:min_s]
        print(f"Using first {min_s} matched samples.")

    if use_log:
        y = np.log1p(y)

    rhos, pvals = [], []

    # compute spearman corrleation between each k-mer frequency and chronological age across all samples
    for col in X.columns:
        rho, p = spearmanr(X[col].values, y)

        if np.isnan(rho):
            rho = 0.0
        if np.isnan(p):
            p = 1.0

        rhos.append(rho)
        pvals.append(p)

    rhos = np.array(rhos)
    pvals = np.array(pvals)

    pvals_safe = np.maximum(pvals, np.nextafter(0, 1))

    # bonferroni correction accounts for multiple hypothesis testing, across all evaluated k-mers
    bonferroni = 0.05 / len(pvals_safe)
    nominal_mask = pvals_safe < 0.05
    bonf_mask = pvals_safe < bonferroni

    n_nominal = int(nominal_mask.sum())
    n_bonf = int(bonf_mask.sum())

    kmer_df = pd.DataFrame({
        "kmer": X.columns,
        "spearman_rho": rhos,
        "p": pvals,
        "neg_log10_p": -np.log10(pvals_safe),
        "nominal_p_lt_0.05": nominal_mask,
        "bonferroni_significant": bonf_mask
    })

    kmer_df = kmer_df.sort_values("spearman_rho", ascending=False)

    all_path = os.path.join(out_dir, f"all_kmer_correlations_{tag}.csv")
    kmer_df.to_csv(all_path, index=False)

    # select k-mers with strongest absolute corrections with age
    top_idx = np.argsort(np.abs(rhos))[::-1][:top_n]

    top_df = pd.DataFrame({
        "kmer": X.columns[top_idx],
        "spearman_rho": rhos[top_idx],
        "p": pvals[top_idx],
        "neg_log10_p": -np.log10(pvals_safe[top_idx]),
        "nominal_p_lt_0.05": nominal_mask[top_idx],
        "bonferroni_significant": bonf_mask[top_idx]
    }).sort_values("spearman_rho", key=lambda s: np.abs(s), ascending=False)

    top_path = os.path.join(out_dir, f"top{top_n}_kmer_correlations_{tag}.csv")
    top_df.to_csv(top_path, index=False)

    print(f"Saved: {all_path}")
    print(f"Saved: {top_path}")
    print(f"Nominal p < 0.05: {n_nominal:,}")
    print(f"Bonferroni significant: {n_bonf:,}")

    fig, axes = plt.subplots(1, 2, figsize=(5.2, 2.6))
    fig.subplots_adjust(wspace=0.42, left=0.10, right=0.98, top=0.88, bottom=0.24)

    n_bins = min(60, max(25, len(rhos) // 50))

    # Panel A
    ax1 = axes[0]
    ax1.hist(
        rhos,
        bins=n_bins,
        color=PAL["blue"],
        edgecolor="white",
        linewidth=0.25,
        alpha=0.88
    )

    ax1.axvline(0, color=PAL["gray"], lw=0.8, ls="--")
    ax1.set_xlabel(r"Spearman $\rho$")
    ax1.set_ylabel("Number of k-mers")
    ax1.set_title(r"A  Correlation with age", loc="left", fontweight="bold")
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))

    ax1.text(
        0.97,
        0.95,
        f"n = {len(rhos):,}",
        transform=ax1.transAxes,
        ha="right",
        va="top",
        fontsize=6.5,
        color=PAL["gray"]
    )

    # Panel B
    ax2 = axes[1]
    ax2.hist(
        pvals_safe,
        bins=n_bins,
        color=PAL["orange"],
        edgecolor="white",
        linewidth=0.25,
        alpha=0.88
    )

    ax2.set_yscale("log")

    ax2.axvline(
        0.05,
        color=PAL["p05col"],
        lw=0.9,
        ls="--",
        label=f"$p < 0.05$: {n_nominal:,}"
    )

    ax2.axvline(
        bonferroni,
        color=PAL["bonfcol"],
        lw=0.9,
        ls=":",
        label=f"Bonferroni: {n_bonf:,}"
    )

    ax2.set_xlabel("$p$-value")
    ax2.set_ylabel("Number of k-mers")
    ax2.set_title("B  $p$-value distribution", loc="left", fontweight="bold")

    ax2.legend(
        frameon=False,
        handlelength=1.4,
        loc="upper right",
        fontsize=6.2
    )

    for ext in ("pdf", "png", "svg"):
        path = os.path.join(out_dir, f"kmer_analysis_2panel_{tag}.{ext}")
        plt.savefig(path)
        print(f"Saved: {path}")

    plt.close()

    return rhos, pvals, kmer_df


# Main prediction function
def run_age_prediction(matrix_file, phenotype_file, tag="run", out_dir="results"):
    """Run reproducible dog age prediction using Ridge regression with true LOOCV."""
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n=== Ridge Age Prediction Model ({tag}) ===")

    # load k-mer feature matrix generated from mkmc
    try:
        X = pd.read_csv(matrix_file, sep="\t", index_col=0)
    except Exception:
        X = pd.read_csv(matrix_file, sep=",", index_col=0)

    print(f"Raw matrix shape: {X.shape}")

    if X.shape[0] > X.shape[1]:
        print("Transposing matrix so samples = rows and k-mers = columns")
        X = X.T

    print(f"Processed: {X.shape[0]} samples × {X.shape[1]} features")

    y = pd.read_csv(phenotype_file, header=None).iloc[:, 0].values

    print(f"Loaded {len(y)} ages")

    if len(y) != X.shape[0]:
        print(
            f"WARNING: sample mismatch. Matrix has {X.shape[0]} samples, "
            f"phenotype has {len(y)} ages."
        )

        min_s = min(len(y), X.shape[0])
        X = X.iloc[:min_s]
        y = y[:min_s]

        print(f"Using first {min_s} matched samples.")

    # optional log transformation to reduce skew in age distribution
    if USE_LOG_TRANSFORM:
        y = np.log1p(y)
        print("Using log1p-transformed ages.")

    if len(y) >= 1000:
        LOO_PRINT_EVERY = 10
    else:
        LOO_PRINT_EVERY = 5

    print("\nRidge settings:")
    print(f"  alpha grid = {RIDGE_ALPHAS[0]:.1e} to {RIDGE_ALPHAS[-1]:.1e}")
    print("  StandardScaler = False")
    print("  alpha selected by RidgeCV inside each LOOCV training fold")

    print(f"\n=== LOOCV: {len(y)} iterations ===")

    # Leave-One-Out Cross Validation,  each sample is used once as the test sample while all, remaining samples are used for training.
    loo = LeaveOneOut()
    y_true_all, y_pred_all = [], []
    alpha_all = []

    for fold, (train_idx, test_idx) in enumerate(loo.split(X)):
        if (fold + 1) % LOO_PRINT_EVERY == 0 or fold == 0:
            print(f"  Iteration {fold + 1}/{len(y)}...")

        # ridge regression reduces overfit in high-dimensional, k-mer feature space through L2 regularization
        model = RidgeCV(alphas=RIDGE_ALPHAS)

        model.fit(X.iloc[train_idx], y[train_idx])
        pred = model.predict(X.iloc[test_idx])[0]

        y_pred_all.append(pred)
        y_true_all.append(y[test_idx][0])
        alpha_all.append(model.alpha_)

    y_true_all = np.array(y_true_all)
    y_pred_all = np.array(y_pred_all)
    alpha_all = np.array(alpha_all)

    r2 = r2_score(y_true_all, y_pred_all)
    rmse = np.sqrt(mean_squared_error(y_true_all, y_pred_all))
    mae = np.mean(np.abs(y_true_all - y_pred_all))
    rho, pval = spearmanr(y_true_all, y_pred_all)

    print("\n=== Overall LOOCV Ridge Results ===")
    print(f"R² = {r2:.3f}")
    print(f"RMSE = {rmse:.3f}")
    print(f"MAE = {mae:.3f}")
    print(f"Spearman rho = {rho:.3f}")
    print(f"p-value = {pval:.3e}")
    print(f"Median selected alpha = {np.median(alpha_all):.4g}")

    # save predictions, residuals, and selected alpha values
    pred_df = pd.DataFrame({
        "y_true": y_true_all,
        "y_pred": y_pred_all,
        "residual": y_pred_all - y_true_all,
        "selected_alpha": alpha_all
    })

    pred_path = os.path.join(out_dir, f"predictions_{tag}.csv")
    pred_df.to_csv(pred_path, index=False)

    print(f"Saved: {pred_path}")

    plot_scatter(y_true_all, y_pred_all, tag=tag, out_dir=out_dir)

    return y_true_all, y_pred_all

# Run
if __name__ == "__main__":
    out_dir = "results"

    matrix_file = "1k_size_6_raw_cs_matrix_norm" # input file
    phenotype_file = "1k_phenotype.txt"  # phenotype file

    y_true, y_pred = run_age_prediction(
        matrix_file,
        phenotype_file,
        tag="blood_ridge_loocv",
        out_dir=out_dir
    )

    rhos, pvals, kmer_df = plot_kmer_correlations(
        matrix_file,
        phenotype_file,
        use_log=False,
        top_n=20,
        out_dir=out_dir,
        tag="blood"
    )
