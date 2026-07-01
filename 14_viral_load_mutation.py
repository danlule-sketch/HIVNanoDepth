#!/usr/bin/env python3

'''

This script analyzes the relationship between plasma viral load and HIV mutation 
burden derived from Nanopore sequencing data. It integrates viral load measurements 
with per-sample mutation summaries, computes total mutation counts across multiple 
drug-resistance categories, performs correlation and regression analyses, and 
generates a log-scaled scatter plot with a fitted regression line.

'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from sklearn.linear_model import LinearRegression

# -----------------------------
# 1. Load data
# -----------------------------
viral_load = pd.read_csv("viral_load.csv")
mutations = pd.read_csv("sequenceSummaries_nanopore.csv")

# Standardise barcode column names
viral_load.rename(columns={"IDNO": "Barcode"}, inplace=True)
mutations.rename(columns={"Sequence Name": "Barcode"}, inplace=True)

# -----------------------------
# 2. Define mutation columns
# -----------------------------
mutation_cols = [
    "PI Major", "PI Accessory", "PR Other",
    "NRTI", "NNRTI", "RT Other",
    "INSTI Major", "INSTI Accessory", "IN Other"
]

# -----------------------------
# 3. Count mutations per row
# -----------------------------
def count_mutations(row):
    total = 0
    for col in mutation_cols:
        val = row.get(col)
        if pd.notna(val) and val not in ["None", "NA", ""]:
            total += len([m for m in val.split(",") if m.strip()])
    return total

mutations["Mutation_Count"] = mutations.apply(count_mutations, axis=1)

# -----------------------------
# 4. Merge with viral load
# -----------------------------
df = pd.merge(
    viral_load,
    mutations[["Barcode", "Mutation_Count"]],
    on="Barcode",
    how="inner"
)

# Drop zeros / invalid values if any
df = df[df["Plasma_Viral _Load(copies/mL)"] > 0]

# -----------------------------
# 5. Statistics
# -----------------------------
spearman_r, spearman_p = spearmanr(
    df["Plasma_Viral _Load(copies/mL)"],
    df["Mutation_Count"]
)

pearson_r, pearson_p = pearsonr(
    np.log10(df["Plasma_Viral _Load(copies/mL)"]),
    df["Mutation_Count"]
)

# -----------------------------
# 6. Regression line
# -----------------------------
X = np.log10(df["Plasma_Viral _Load(copies/mL)"]).values.reshape(-1, 1)
y = df["Mutation_Count"].values

reg = LinearRegression()
reg.fit(X, y)
y_pred = reg.predict(X)

# -----------------------------
# 7. Plot
# -----------------------------
plt.figure(figsize=(8, 6))

sns.set_style("whitegrid")
sns.set_context("notebook", font_scale=1.3)

# Scatter
sns.scatterplot(
    data=df,
    x="Plasma_Viral _Load(copies/mL)",
    y="Mutation_Count",
    s=100,
    color="steelblue"
)

# Regression line
plt.plot(
    10**X.flatten(),  # back-transform log10 for plotting
    y_pred,
    color="darkred",
    linewidth=2,
    label="Linear Regression"
)

plt.xscale("log")

plt.xlabel("Plasma Viral Load (copies/mL, $\\log_{10}$)", fontsize=20, fontweight="bold")
plt.ylabel("Total Mutation Count", fontsize=20, fontweight="bold")
plt.title("Viral Load vs Mutation Burden", fontsize=25, fontweight="bold")

# Annotate statistics
stats_text = (
    f"Spearman ρ = {spearman_r:.2f} (p = {spearman_p:.3g})\n"
    f"Pearson r = {pearson_r:.2f} (p = {pearson_p:.3g})"
)

plt.text(
    0.05, 0.95,
    stats_text,
    transform=plt.gca().transAxes,
    fontsize=12,
    verticalalignment="top",
    fontweight="bold",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)

plt.legend(fontsize=12, frameon=True)
plt.tight_layout()
plt.savefig("viral_load_vs_mutation_count_with_regression.png", dpi=300)
plt.show()
