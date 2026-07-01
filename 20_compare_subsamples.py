#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

'''
This script visualizes how variant-calling performance (F1 score) changes across 
different BAM subsampling proportions i.e. 0.005, 0.01, 0.02, 0.05, 0.1

It produces a log-scaled boxplot of F1 scores at each sampling proportion, with:

	Log-aware box widths
	Overlaid individual data points
	Sample size annotations
	Publication-quality styling

This creates a box plot designed to assess whether variant-calling accuracy is robust 
to sequencing depth reduction.

'''

# -----------------------------
# 1. Load data
# -----------------------------
df = pd.read_csv("comparison_results/metrics_happy.tsv", sep="\t")

df.columns = [c.lower() for c in df.columns]

required_cols = ["proportion", "f1", "barcode"]
for c in required_cols:
    if c not in df.columns:
        raise KeyError(f"Expected column '{c}' not found")

df["proportion"] = df["proportion"].astype(float)
df = df.sort_values("proportion")

# -----------------------------
# 2. Prepare grouped data
# -----------------------------
grouped = df.groupby("proportion")["f1"]

proportions = np.array(sorted(grouped.groups.keys()))
data = [grouped.get_group(p).values for p in proportions]
counts = grouped.size()

# -----------------------------
# 3. Compute log-aware box widths
# -----------------------------
log_props = np.log10(proportions)
log_diffs = np.diff(log_props)
log_diffs = np.append(log_diffs, log_diffs[-1])
widths = (10 ** (log_props + log_diffs * 0.3)) - proportions

# -----------------------------
# 4. Plot setup (journal style)
# -----------------------------
sns.set_style("white")
plt.figure(figsize=(14, 8))

plt.rcParams.update({
    "axes.labelsize": 20,
    "axes.titlesize": 24,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18
})

# -----------------------------
# 5. Numeric log-x boxplots
# -----------------------------
bp = plt.boxplot(
    data,
    positions=proportions,
    widths=widths,
    patch_artist=True,
    showfliers=False
)

# Style boxes
for box in bp["boxes"]:
    box.set(facecolor="#9ecae1", edgecolor="black", linewidth=1.5)

for element in ["medians", "whiskers", "caps"]:
    for line in bp[element]:
        line.set(color="black", linewidth=1.5)

# Optional: overlay raw points (very light)
for p in proportions:
    y = grouped.get_group(p)
    x = np.random.normal(p, p * 0.015, size=len(y))
    plt.scatter(x, y, s=18, color="black", alpha=0.4, zorder=3)

# -----------------------------
# 6. Annotate sample counts
# -----------------------------
ymin, ymax = plt.ylim()
y_text = ymin + 0.02 * (ymax - ymin)

for p in proportions:
    plt.text(
        p,
        y_text,
        f"n={counts[p]}",
        ha="center",
        va="bottom",
        fontsize=14
    )

# -----------------------------
# 7. Formatting
# -----------------------------
plt.xscale("log")
plt.xticks(
    [0.005, 0.01, 0.02, 0.05, 0.1],
    ["0.005", "0.01", "0.02", "0.05", "0.1"]
)

plt.xlabel("Sampling proportion", fontweight="bold")
plt.ylabel("F1 score", fontweight="bold")
plt.title(
    "Variant Calling Performance Is Invariant to BAM Subsampling Depth",
    fontweight="bold"
)

# Clean spines (Nature/Bioinformatics style)
for spine in ["top", "right"]:
    plt.gca().spines[spine].set_visible(False)

plt.tight_layout()

# -----------------------------
# 8. Save
# -----------------------------
plt.savefig(
    "comparison_results/F1_distribution_numeric_logx.png",
    dpi=600
)
plt.savefig(
    "comparison_results/F1_distribution_numeric_logx.pdf"
)

plt.show()
