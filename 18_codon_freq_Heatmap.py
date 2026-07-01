#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re

'''

This script generates a heatmap visualization of Nanopore variant frequencies 
at codon positions that were also detected and confirmed by Sanger sequencing.

The script inputs a csv file "global_comparison_summary.csv" generated from the
output of the script "17_codon_bias_from_bam" which previous input folders
/fasta_pr_references_fowardclear and /minimap2_results2_medaka. The final output is a 
PNG heatmap showing Nanopore frequency (%) per barcode and codon position.

'''

# -------------------------------
# LOAD DATA
# -------------------------------
df = pd.read_csv("global_comparison_summary.csv")
subtype_df = pd.read_csv("sequenceSummaries.csv", on_bad_lines="skip")

# -------------------------------
# CLEAN BARCODE FUNCTION
# -------------------------------
def clean_bc(x):
    if pd.isna(x):
        return None
    x = str(x).lower()

    # extract core barcode (first meaningful token)
    m = re.search(r"(barcode\d+|[a-z]*\d+)", x)
    if not m:
        return None

    return m.group(1)

# apply cleaning
df["bc"] = df["barcode"].apply(clean_bc)
subtype_df["bc"] = subtype_df["Sequence Name"].apply(clean_bc)

# -------------------------------
# CLEAN SUBTYPE
# -------------------------------
subtype_df["subtype"] = (
    subtype_df["Subtype (%)"]
    .astype(str)
    .str.extract(r"([A-Za-z0-9]+)")
)

subtype_df = subtype_df.dropna(subset=["bc", "subtype"])
subtype_df = subtype_df[["bc", "subtype"]]

# -------------------------------
# FILTER DATA
# -------------------------------
sanger_df = df[df["Detected_by_Sanger"] == True].copy()

match_df = sanger_df[
    sanger_df["Sanger_codon"] == sanger_df["Nanopore_codon"]
].copy()

# -------------------------------
# MERGE SUBTYPE (NOW ROBUST)
# -------------------------------
match_df = match_df.merge(subtype_df, on="bc", how="left")

# debug check (IMPORTANT)
print("Subtype match rate:", match_df["subtype"].notna().mean())

match_df["subtype"] = match_df["subtype"].fillna("UNK")

# -------------------------------
# PIVOT
# -------------------------------
heat_df = match_df.pivot_table(
    index="barcode",
    columns="position",
    values="nano_freq(%)",
    aggfunc="max",
    fill_value=0
)

# -------------------------------
# ORDER
# -------------------------------
order_df = (
    match_df[["barcode", "subtype"]]
    .drop_duplicates()
    .sort_values(["subtype", "barcode"])
)

heat_df = heat_df.loc[order_df["barcode"]]
subtype_series = order_df.set_index("barcode")["subtype"]

# -------------------------------
# FILTER OUT EMPTY SUBTYPES
# -------------------------------
valid = subtype_series != "UNK"
subtype_series = subtype_series[valid]

# align heatmap too
heat_df = heat_df.loc[subtype_series.index]

# -------------------------------
# COLORS
# -------------------------------
unique_subtypes = subtype_series.unique()
palette = sns.color_palette("tab20", len(unique_subtypes))
color_map = dict(zip(unique_subtypes, palette))

# -------------------------------
# PLOT
# -------------------------------
plt.figure(figsize=(26, 16))
ax = sns.heatmap(
    heat_df,
    cmap="YlGnBu",
    cbar_kws={"label": "Nanopore Frequency (%)"}
)

# -------------------------------
# AXES CLEANUP (IMPORTANT)
# -------------------------------
ax.set_xlabel("Codon Position", fontsize=22, fontweight="bold")
ax.set_ylabel("Barcode", fontsize=22, fontweight="bold")
ax.set_title(
    "Nanopore Frequencies at Sanger Codon Positions (Grouped by Subtype)",
    fontsize=26,
    fontweight="bold"
)

# reduce x-axis clutter
n_x = len(heat_df.columns)
step = max(1, n_x // 10)

ax.set_xticks(np.arange(n_x)[::step])
ax.set_xticklabels(
    heat_df.columns[::step],
    rotation=90,
    fontsize=12,
    fontweight="bold"
)

# y-axis reduction
n_y = len(heat_df.index)
y_step = max(1, n_y // 30)

ax.set_yticks(np.arange(n_y)[::y_step])
ax.set_yticklabels(
    heat_df.index[::y_step],
    fontsize=10,
    fontweight="bold"
)

# -------------------------------
# SUBTYPE BRACKETS (FINAL FIXED)
# -------------------------------
barcode_to_y = {b: i for i, b in enumerate(heat_df.index)}

group_df = subtype_series.reset_index()
group_df["y"] = group_df["barcode"].map(barcode_to_y)
group_df = group_df.dropna().sort_values("y")

groups = group_df.groupby("subtype")

x_bracket = heat_df.shape[1] + 1.5
x_text = heat_df.shape[1] + 4.5
lw = 5

for subtype, g in groups:
    y = g["y"].values
    y0, y1 = int(y.min()), int(y.max()) + 1
    yc = (y0 + y1) / 2

    color = color_map[subtype]

    ax.plot([x_bracket, x_bracket], [y0, y1], color=color, lw=lw, clip_on=False)
    ax.plot([x_bracket - 2, x_bracket + 2], [y0, y0], color=color, lw=lw, clip_on=False)
    ax.plot([x_bracket - 2, x_bracket + 2], [y1, y1], color=color, lw=lw, clip_on=False)

    ax.text(
        x_text, yc, subtype,
        fontsize=18,
        fontweight="bold",
        color=color,
        va="center",
        clip_on=False
    )

# -------------------------------
# COLORBAR
# -------------------------------
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=14)
cbar.set_label("Nanopore Frequency (%)", fontsize=16, fontweight="bold")

plt.tight_layout()
plt.savefig("FIXED_subtype_heatmap.png", dpi=300, bbox_inches="tight")
plt.show()