#!/usr/bin/env python3

'''
This script integrates sequencing depth data and plasma viral load measurements to 
evaluate the relationship between viral load and sequencing coverage. It computes 
summary depth statistics per sample, merges them with viral load values based on 
matched sample identifiers, performs correlation and regression analyses, and 
generates a log–log scatter plot with regression line and statistical annotations.

The output figure includes both Pearson and Spearman correlation metrics.

'''

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy import stats
import numpy as np
import re

# -------------------------------
# Paths
# -------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

DEPTH_FILE = os.path.join(BASE_DIR, "all_depths.txt")
VIRAL_LOAD_FILE = os.path.join(BASE_DIR, "viral_load.csv")

# -------------------------------
# Load depth data
# -------------------------------
depths = pd.read_csv(
    DEPTH_FILE,
    sep=r"\s+",
    header=None,
    names=["BarcodeFile", "Sample", "Position", "Depth"]
)

# Keep only positive depths
depths = depths[depths["Depth"] > 0]

# Summarise depth per sample
depth_summary = (
    depths
    .groupby("Sample")["Depth"]
    .agg(
        mean_depth="mean",
        median_depth="median",
        min_depth="min",
        max_depth="max"
    )
    .reset_index()
)

# -------------------------------
# Load viral load data
# -------------------------------
viral_load = pd.read_csv(VIRAL_LOAD_FILE)

# Extract barcode number (01, 02, ...)
viral_load["barcode_num"] = (
    viral_load["IDNO"]
    .str.extract(r"(\d+)")
    .astype(int)
)

# Match depth sample format: V24_05_010
viral_load["Sample"] = viral_load["barcode_num"].apply(
    lambda x: f"V24_05_{x:03d}"
)

viral_load = viral_load.rename(
    columns={"Plasma_Viral _Load(copies/mL)": "Viral_Load"}
)

# Keep only positive viral loads
viral_load = viral_load[viral_load["Viral_Load"] > 0]

# -------------------------------
# Merge depth + viral load
# -------------------------------
merged = pd.merge(
    viral_load,
    depth_summary,
    on="Sample",
    how="inner"
)

print(f"Merged samples: {len(merged)}")

if len(merged) < 2:
    raise ValueError(
        "Not enough matched samples for correlation analysis.\n"
        "Check that sample IDs match between viral_load.csv and all_depths.txt."
    )

# -------------------------------
# Statistics
# -------------------------------
x = merged["Viral_Load"]
y = merged["mean_depth"]

pearson_r, pearson_p = stats.pearsonr(x, y)
spearman_r, spearman_p = stats.spearmanr(x, y)

slope, intercept, r_value, p_value, std_err = stats.linregress(
    np.log10(x),
    np.log10(y)
)

# -------------------------------
# Plot
# -------------------------------
plt.figure(figsize=(9, 7))

ax = sns.regplot(
    data=merged,
    x="Viral_Load",
    y="mean_depth",
    scatter_kws={"s": 80, "alpha": 0.85},
    line_kws={"color": "red", "linewidth": 2}
)

# Explicit log scales
ax.set_xscale("log")
ax.set_yscale("log")

# Axis labels (larger + bold)
ax.set_xlabel(
    "Plasma Viral Load (copies/mL)",
    fontsize=20,
    fontweight="bold"
)
ax.set_ylabel(
    "Mean Sequencing Depth",
    fontsize=20,
    fontweight="bold"
)

# Title (larger + bold)
ax.set_title(
    "Viral Load vs Sequencing Depth",
    fontsize=25,
    fontweight="bold",
    pad=12
)

# Increase tick label size
ax.tick_params(axis="both", which="major", labelsize=16)
ax.tick_params(axis="both", which="minor", labelsize=16)

# Stats annotation
stats_text = (
    f"n = {len(merged)}\n"
    f"Pearson r = {pearson_r:.2f} (p = {pearson_p:.2e})\n"
    f"Spearman ρ = {spearman_r:.2f} (p = {spearman_p:.2e})\n"
    f"Log–log R² = {r_value**2:.2f}"
)

ax.text(
    0.05, 0.95,
    stats_text,
    transform=ax.transAxes,
    ha="left",
    va="top",
    fontsize=16,
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.85)
)

plt.tight_layout()

# -------------------------------
# Save
# -------------------------------
output_file = os.path.join(
    BASE_DIR,
    "viral_load_vs_mean_depth.png"
)
plt.savefig(output_file, dpi=300)
plt.show()
