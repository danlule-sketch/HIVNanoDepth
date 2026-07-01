#!/usr/bin/env python


'''
This script reads per-base sequencing depth data, computes the average sequencing 
depth per sample, applies outlier mitigation using IQR-based winsorization, and 
generates a scatter plot showing average depth across samples 
on an absolute (linear) scale.

The visualization is optimized for large numbers of samples by skipping crowded x-axis 
labels and applying consistent styling suitable for manuscripts.

The script input file is all_depths.txt (change line 32 appropriately to point to your path)
the script creates a plot average_depth_all_samples_absolute_scale_skipped_labels.png
The folder 

'''
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import MaxNLocator, AutoMinorLocator

# -----------------------------
# Style for publication
# -----------------------------
sns.set_theme(
    style="ticks",
    context="paper",
    font_scale=1.2
)

# -----------------------------
# Read depth file
# -----------------------------
depths = pd.read_csv(
    "depths_minimap2/all_depths.txt", # Change the path to your input depth file location
    sep=" ",
    header=None,
    names=["Sample", "Chromosome", "Position", "Depth"]
)

# -----------------------------
# Compute average depth per sample
# -----------------------------
summary = (
    depths
    .groupby("Sample", as_index=False)
    .agg(AvgDepth=("Depth", "mean"))
)

# -----------------------------
# Outlier handling (IQR winsorization)
# -----------------------------
def winsorize(series):
    q1, q3 = series.quantile([0.25, 0.75])
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    return series.clip(lower, upper)

summary["Depth_w"] = winsorize(summary["AvgDepth"])

# -----------------------------
# Clean and sort sample labels
# -----------------------------
summary["SampleLabel"] = summary["Sample"].str.split(".").str[0]
summary = summary.sort_values("SampleLabel").reset_index(drop=True)
summary["SampleID"] = range(1, len(summary) + 1)

# -----------------------------
# Figure size
# -----------------------------
fig, ax = plt.subplots(figsize=(9.5, 6.5))

# -----------------------------
# Depth plot (axes swapped)
# -----------------------------
sns.scatterplot(
    data=summary,
    x="SampleID",
    y="Depth_w",
    ax=ax,
    color="#C44E52",
    s=40
)

# -----------------------------
# Labels and title (ABSOLUTE SCALE)
# -----------------------------
ax.set_xlabel("Sample", fontsize=22)
ax.set_ylabel("Average Depth", fontsize=22)

ax.set_title(
    "Average Depth per Sample",
    fontsize=26,
    fontweight="bold",
    pad=15
)

# -----------------------------
# Skip crowded x-axis labels
# -----------------------------
N = 5  # show every N-th sample

xticks = summary["SampleID"]
xlabels = summary["SampleLabel"].where(summary.index % N == 0, "")

ax.set_xticks(xticks)
ax.set_xticklabels(
    xlabels,
    rotation=90,
    ha="center",
    fontsize=14
)

# -----------------------------
# Y-axis tick control (linear scale)
# -----------------------------
ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.tick_params(axis="y", which="major", labelsize=16, length=6)
ax.tick_params(axis="y", which="minor", length=3)

# -----------------------------
# Clean look
# -----------------------------
sns.despine(trim=True)

# -----------------------------
# Layout
# -----------------------------
plt.subplots_adjust(
    left=0.12,
    right=0.97,
    bottom=0.35,
    top=0.90
)

# -----------------------------
# Save figure
# -----------------------------
plt.savefig(
    "average_depth_all_samples_absolute_scale_skipped_labels.png",
    dpi=300
)
plt.show()
