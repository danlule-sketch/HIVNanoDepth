#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix, cohen_kappa_score
import re


'''
This script evaluates HIV subtype concordance between Nanopore and Sanger sequencing 
results using per-sample subtype calls.

The script inputs 2 csv files generated from both the Sanger and Nanopore HIVDR
output from the Stanford HIVDR database; "sequenceSummaries_nanopore.csv" and
"sequenceSummaries_sanger.csv"

The output is a confusion matrix figure

'''
# =========================
# 1. Load data
# =========================
nano = pd.read_csv("sequenceSummaries_nanopore.csv")
sanger = pd.read_csv("sequenceSummaries_sanger.csv")

# =========================
# 2. Harmonize sample IDs
# =========================
nano["sample_id"] = nano["Sequence Name"].str.strip()
sanger["sample_id"] = sanger["Sequence Name"].str.extract(r"(Barcode\d+)")

# =========================
# 3. Normalize subtype strings
# =========================
def normalize_subtype(x):
    if pd.isna(x):
        return "Unknown"
    m = re.match(r"([A-Z]\d?)", x.strip())
    return m.group(1) if m else "Other"

nano["Subtype_norm"] = nano["Subtype (%)"].apply(normalize_subtype)
sanger["Subtype_norm"] = sanger["Subtype (%)"].apply(normalize_subtype)

# =========================
# 4. Merge datasets
# =========================
merged = nano.merge(
    sanger,
    on="sample_id",
    suffixes=("_nano", "_sanger")
)

print(f"Matched samples: {merged.shape[0]}")

# =========================
# 5. Subtype concordance stats
# =========================
subtype_df = merged[
    ["sample_id", "Subtype_norm_nano", "Subtype_norm_sanger"]
]

percent_agreement = (
    subtype_df["Subtype_norm_nano"] ==
    subtype_df["Subtype_norm_sanger"]
).mean() * 100

kappa = cohen_kappa_score(
    subtype_df["Subtype_norm_nano"],
    subtype_df["Subtype_norm_sanger"]
)

print("\nSubtype concordance:")
print(f"Percent agreement: {percent_agreement:.2f}%")
print(f"Cohen's kappa: {kappa:.2f}")

# =========================
# 6. Save subtype-discordant samples
# =========================
discordant = subtype_df[
    subtype_df["Subtype_norm_nano"] != subtype_df["Subtype_norm_sanger"]
]

discordant.to_csv(
    "subtype_discordant_samples.csv",
    index=False
)

print(
    f"Saved {discordant.shape[0]} subtype-discordant samples "
    f"to subtype_discordant_samples.csv"
)

# =========================
# 7. Confusion matrix (bold, large font) with Cohen's kappa
# =========================
labels = sorted(
    set(subtype_df["Subtype_norm_nano"]) |
    set(subtype_df["Subtype_norm_sanger"])
)

# Compute confusion matrix
cm = confusion_matrix(
    subtype_df["Subtype_norm_sanger"],
    subtype_df["Subtype_norm_nano"],
    labels=labels
)

plt.figure(figsize=(8, 7))
ax = sns.heatmap(
    cm,
    annot=True,
    fmt="d",
    cmap="Blues",
    cbar=True,
    xticklabels=labels,
    yticklabels=labels,
    annot_kws={"size": 18, "weight": "bold"}
)

# Axis labels
ax.set_xlabel("Nanopore Subtype", fontsize=18, fontweight="bold")
ax.set_ylabel("Sanger Subtype", fontsize=18, fontweight="bold")

# Tick labels
ax.set_xticklabels(
    ax.get_xticklabels(),
    fontsize=16,
    fontweight="bold",
    rotation=0,
    ha="right"
)
ax.set_yticklabels(
    ax.get_yticklabels(),
    fontsize=16,
    fontweight="bold",
    rotation=0
)

# Title
ax.set_title(
    "Subtype Concordance: Nanopore vs Sanger",
    fontsize=20,
    fontweight="bold",
    pad=20
)

# Cohen's kappa + percent agreement displayed below heatmap
kappa_text = f"Cohen’s κ = {kappa:.2f}\nPercent agreement = {percent_agreement:.1f}%"
ax.text(
    0.5,
    -0.25,  # position below heatmap
    kappa_text,
    transform=ax.transAxes,
    ha="center",
    va="top",
    fontsize=16,
    fontweight="bold"
)

plt.tight_layout()
plt.savefig(
    "nanopore_sanger_subtype_confusion.png",
    dpi=300
)
plt.show()
