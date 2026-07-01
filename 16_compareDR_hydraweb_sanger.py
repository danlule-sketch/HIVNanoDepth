#!/usr/bin/env python3

'''
This script performs a position-level comparison of HIV mutation calls derived from 
Sanger sequencing and Oxford Nanopore sequencing. It quantifies how consistently 
mutation positions are detected between the two platforms across matched samples 
and summarizes agreement using descriptive statistics, hypothesis testing, and 
publication-quality visualizations.

The analysis operates at the (gene, position) level and does not require identical 
mutation identities—only positional concordance.

The script expects the following directory structure:

fastq_filtered/
├── mutation-list_Sanger_stanford/
│   ├── MutationList_sample1.csv
│   ├── MutationList_sample2.csv
│   └── ...
├── mutation-list_nanopore_stanford/
│   ├── MutationList_sample1.csv
│   ├── MutationList_sample2.csv
│   └── ...

'''

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

# ==========================
# PATHS (READY TO RUN)
# ==========================
BASE_DIR = "/Users/daniellulebugembe/Documents/PhD/reading7_simulation/minimap2_results2_medaka/fastq_filtered"

SANGER_DIR = os.path.join(BASE_DIR, "mutation-list_Sanger_stanford")
NANOPORE_DIR = os.path.join(BASE_DIR, "mutation-list_nanopore_stanford")

# ==========================
# RESULTS DIRECTORY
# ==========================
RESULTS_DIR = os.path.join(BASE_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

OUT_TABLE = os.path.join(RESULTS_DIR, "position_comparison_table.csv")
OUT_STATS = os.path.join(RESULTS_DIR, "position_concordance_stats.csv")

OUT_SCATTER = os.path.join(RESULTS_DIR, "position_concordance_scatter.png")
OUT_BOXPLOT_COUNTS = os.path.join(RESULTS_DIR, "position_counts_boxplot.png")
OUT_VIOLIN_JACCARD = os.path.join(RESULTS_DIR, "jaccard_summary_violin.png")

# ==========================
# JOURNAL-STYLE PLOT HELPER
# ==========================
def journal_style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.tick_params(width=1.2, labelsize=14)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")

# ==========================
# LOAD MUTATION LIST FILES
# ==========================
def load_mutation_folder(folder, source):
    records = []

    for file in glob.glob(os.path.join(folder, "*.csv")):
        df = pd.read_csv(file)

        if "Sequence Name" in df.columns:
            sample = df["Sequence Name"].iloc[0].split("_")[0]
        else:
            sample = os.path.basename(file)\
                .replace("MutationList_", "")\
                .replace(".csv", "")

        df = df[["Gene", "Position"]].drop_duplicates()
        df["Sample"] = sample
        df["Source"] = source
        records.append(df)

    return pd.concat(records, ignore_index=True)

print("Loading Sanger mutation lists...")
sanger_df = load_mutation_folder(SANGER_DIR, "Sanger")

print("Loading Nanopore mutation lists...")
nano_df = load_mutation_folder(NANOPORE_DIR, "Nanopore")

# ==========================
# MERGE & DEDUPLICATE
# ==========================
all_positions = pd.concat([sanger_df, nano_df], ignore_index=True)
all_positions = all_positions.drop_duplicates(
    subset=["Sample", "Gene", "Position", "Source"]
)

# ==========================
# POSITION COMPARISON
# ==========================
comparison_rows = []

for sample in sorted(all_positions["Sample"].unique()):
    sanger_pos = set(
        all_positions.query(
            "Sample == @sample and Source == 'Sanger'"
        )[["Gene", "Position"]].apply(tuple, axis=1)
    )

    nano_pos = set(
        all_positions.query(
            "Sample == @sample and Source == 'Nanopore'"
        )[["Gene", "Position"]].apply(tuple, axis=1)
    )

    for pos in sanger_pos | nano_pos:
        if pos in sanger_pos and pos in nano_pos:
            status = "Shared"
        elif pos in sanger_pos:
            status = "Sanger only"
        else:
            status = "Nanopore only"

        comparison_rows.append({
            "Sample": sample,
            "Gene": pos[0],
            "Position": pos[1],
            "Status": status
        })

comparison_df = pd.DataFrame(comparison_rows)
comparison_df.to_csv(OUT_TABLE, index=False)

# ==========================
# PER-SAMPLE STATISTICS
# ==========================
stats = []

for sample in sorted(comparison_df["Sample"].unique()):
    sub = comparison_df.query("Sample == @sample")

    sanger_n = (sub["Status"] != "Nanopore only").sum()
    nano_n = (sub["Status"] != "Sanger only").sum()
    shared_n = (sub["Status"] == "Shared").sum()
    union_n = sub.shape[0]

    jaccard = shared_n / union_n if union_n else 0

    stats.append({
        "Sample": sample,
        "Sanger_positions": sanger_n,
        "Nanopore_positions": nano_n,
        "Shared_positions": shared_n,
        "Jaccard": jaccard
    })

stats_df = pd.DataFrame(stats)
stats_df.to_csv(OUT_STATS, index=False)

# ==========================
# STATISTICS
# ==========================
w_stat, p_value = wilcoxon(
    stats_df["Nanopore_positions"],
    stats_df["Sanger_positions"]
)

mean_j = stats_df["Jaccard"].mean()
std_j = stats_df["Jaccard"].std()
n = len(stats_df)

ci_low = mean_j - 1.96 * std_j / np.sqrt(n)
ci_high = mean_j + 1.96 * std_j / np.sqrt(n)

# ==========================
# SCATTER: COUNTS
# ==========================
plt.figure(figsize=(14, 6))
ax = plt.gca()

ax.scatter(
    stats_df["Sanger_positions"],
    stats_df["Nanopore_positions"],
    s=60,
    color="#1f77b4",
    alpha=0.85
)

max_val = max(
    stats_df["Sanger_positions"].max(),
    stats_df["Nanopore_positions"].max()
)
ax.plot([0, max_val], [0, max_val], "--", color="black", linewidth=1)

ax.set_xlabel("Sanger mutation positions", fontsize=18, fontweight="bold")
ax.set_ylabel("Nanopore mutation positions", fontsize=18, fontweight="bold")

ax.set_title(
    "Position-level Concordance Between Sanger and Nanopore",
    fontsize=20,
    fontweight="bold",
    pad=12
)

ax.text(
    1, 1,
    f"Wilcoxon p = {p_value:.2e}",
    transform=ax.transAxes,
    fontsize=14,
    fontweight="bold",
    va="top"
)

journal_style(ax)
plt.tight_layout()
plt.savefig(OUT_SCATTER, dpi=300)
plt.close()

# ==========================
# BOXPLOT: COUNTS
# ==========================
plt.figure(figsize=(9, 6))
ax = plt.gca()

ax.boxplot(
    [stats_df["Sanger_positions"], stats_df["Nanopore_positions"]],
    widths=0.5,
    patch_artist=True,
    boxprops=dict(facecolor="#d9d9d9", linewidth=1.5),
    medianprops=dict(color="black", linewidth=1.5),
    whiskerprops=dict(linewidth=1.5),
    capprops=dict(linewidth=1.5)
)

ax.set_xticks([1, 2])
ax.set_xticklabels(["Sanger", "Nanopore"], fontsize=20, fontweight="bold")
ax.set_ylabel("Number of mutation positions", fontsize=18, fontweight="bold")

ax.set_title(
    "Comparison of Detected Mutation Positions",
    fontsize=20,
    fontweight="bold",
    pad=12
)

ax.text(
    1, 1,
    f"Wilcoxon p = {p_value:.2e}",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=14,
    fontweight="bold"
)

journal_style(ax)
plt.tight_layout()
plt.savefig(OUT_BOXPLOT_COUNTS, dpi=300)
plt.close()

# ==========================
# VIOLIN: JACCARD SIMILARITY
# ==========================
plt.figure(figsize=(10, 8))
ax = plt.gca()

parts = ax.violinplot(
    stats_df["Jaccard"],
    showmeans=False,
    showmedians=False,
    showextrema=False
)

for pc in parts["bodies"]:
    pc.set_facecolor("#9ecae1")
    pc.set_edgecolor("black")
    pc.set_alpha(0.9)

ax.boxplot(
    stats_df["Jaccard"],
    widths=0.15,
    patch_artist=True,
    boxprops=dict(facecolor="white", linewidth=1.5),
    medianprops=dict(color="black", linewidth=1.5),
    whiskerprops=dict(linewidth=1.2),
    capprops=dict(linewidth=1.2)
)

ax.set_xticks([1])
ax.set_xticklabels(["Sanger vs Nanopore"], fontsize=18, fontweight="bold")
ax.set_ylabel("Jaccard similarity", fontsize=18, fontweight="bold")

ax.set_title(
    "Overall Mutation Concordance Across Samples",
    fontsize=20,
    fontweight="bold",
    pad=12
)

ax.text(
    1, 1,
    f"Mean = {mean_j:.2f}\n95% CI [{ci_low:.2f}, {ci_high:.2f}]",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=14,
    fontweight="bold"
)

journal_style(ax)
plt.tight_layout()
plt.savefig(OUT_VIOLIN_JACCARD, dpi=300)
plt.close()

print("\nOutputs written to:")
print(RESULTS_DIR)
