#!/usr/bin/env python3
import os
import glob
import pysam
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from Bio import SeqIO
from scipy.stats import pearsonr

'''

This script performs a codon-level comparison between Sanger consensus sequences 
and Nanopore sequencing data across multiple barcodes.

For each barcode, it:

	Extracts codons from the Sanger consensus FASTA.
	Reconstructs codon frequencies from Nanopore BAM alignments.
	Compares codons and frequencies between platforms.
	Produces per-barcode and global summary tables and plots.
	Computes overall recovery and correlation statistics.
	
The goal is to quantify how well Nanopore sequencing recovers Sanger-confirmed codons, 
both per sample and globally.

The script requires files from 2 folders;
	1. /fasta_pr_references_fowardclear
	2. /minimap2_results2_medaka

The output is a csv "global_comparison_summary.csv" require as the input for the
script "18_codon_freq_Heatmap.py" and a figure "global_codon_freq_comparison.png"

'''

# -----------------------------
# Directories (edit as needed)
# -----------------------------
SANGER_DIR = "/Users/daniellulebugembe/Documents/PhD/reading7_simulation/fasta_pr_references_fowardclear"
BAM_DIR    = "/Users/daniellulebugembe/Documents/PhD/reading7_simulation/minimap2_results2_medaka"
OUTPUT_DIR = "./sanger_nanopore_comparison_full"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# -----------------------------
# Codon extraction functions
# -----------------------------
def extract_sanger_codons(fasta_file):
    """Extract codons in frame 0 from Sanger consensus FASTA."""
    record = SeqIO.read(fasta_file, "fasta")
    seq = str(record.seq).upper()
    rows = []
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        if all(b in "ACGT" for b in codon):
            rows.append({"position": i+1, "Sanger_codon": codon})
    return pd.DataFrame(rows)

def reconstruct_nanopore_codons_fast(bam_file):
    """Efficiently reconstruct codons from the BAM."""
    bam = pysam.AlignmentFile(bam_file, "rb")
    contig = bam.references[0]
    ref_to_bases = defaultdict(list)

    for read in bam.fetch(contig):
        if read.is_unmapped:
            continue
        seq = read.query_sequence
        if seq is None:
            continue
        ref_pos = read.reference_start
        qpos = 0
        for cig_op, length in read.cigartuples:
            if cig_op == 0:  # match/mismatch
                for i in range(length):
                    if qpos+i < len(seq):
                        base = seq[qpos+i].upper()
                        if base in "ACGT":
                            ref_to_bases[ref_pos+i].append(base)
                ref_pos += length
                qpos += length
            elif cig_op == 1:  # insertion
                qpos += length
            elif cig_op in (2, 3):  # deletion or skip
                ref_pos += length
            else:
                qpos += length

    codon_counts = defaultdict(lambda: defaultdict(int))
    all_positions = sorted(ref_to_bases.keys())
    for pos in all_positions:
        frame_start = pos - (pos % 3)
        triple = []
        for offset in range(3):
            p = frame_start + offset
            if p in ref_to_bases:
                triple.append(ref_to_bases[p])
            else:
                triple = []
                break
        if len(triple)==3:
            for bases in zip(*triple):
                cod = "".join(bases)
                if len(cod)==3:
                    codon_counts[frame_start][cod]+=1

    rows=[]
    for p,d in codon_counts.items():
        total = sum(d.values())
        for cod,count in d.items():
            freq = count/total*100
            rows.append({"position":p+1,"Nanopore_codon":cod,
                         "nano_count":count,"nano_freq(%)":round(freq,2)})
    return pd.DataFrame(rows)

# -----------------------------
# Comparison and output
# -----------------------------
def compare_and_save(barcode, sanger_df, nano_df):
    merged = pd.merge(sanger_df, nano_df,
                      left_on=["position","Sanger_codon"],
                      right_on=["position","Nanopore_codon"],
                      how="outer")

    merged["Sanger_codon"]   = merged["Sanger_codon"].fillna("-")
    merged["Nanopore_codon"] = merged["Nanopore_codon"].fillna("-")
    merged["nano_freq(%)"]   = merged["nano_freq(%)"].fillna(0)
    merged["Detected_by_Sanger"]   = merged["Sanger_codon"] != "-"
    merged["Detected_by_Nanopore"] = merged["Nanopore_codon"] != "-"

    # Add Sanger frequency as 100% for consensus
    merged["frequency(%)_sanger"] = merged["Detected_by_Sanger"].astype(int) * 100

    out_csv = os.path.join(OUTPUT_DIR, f"{barcode}_comparison.csv")
    merged.to_csv(out_csv,index=False)
    return merged

def plot_comparison(barcode, merged):
    plt.figure(figsize=(20,6))
    labels = merged["position"].astype(str) + "_" + \
             merged["Sanger_codon"] + "/" + merged["Nanopore_codon"]
    x = range(len(labels))
    plt.bar(x, merged["nano_freq(%)"], color="orange")
    plt.xticks(x, labels, rotation=90, fontsize=6)
    plt.ylabel("Nanopore freq (%)")
    plt.title(f"{barcode} Nanopore frequencies (Sanger vs Nano)")
    out_png = os.path.join(OUTPUT_DIR, f"{barcode}_comparison.png")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

# -----------------------------
# Main per-barcode loop
# -----------------------------
all_merged = []

sanger_fastas = glob.glob(os.path.join(SANGER_DIR, "*.fasta"))
barcodes = sorted({os.path.basename(x).split("_")[0] for x in sanger_fastas})

for bc in barcodes:
    print(f"Processing barcode {bc}...")
    sanger_file = glob.glob(os.path.join(SANGER_DIR, f"{bc}*.fasta"))
    bam_file    = os.path.join(BAM_DIR, f"{bc}.sorted.bam")
    if not sanger_file:
        print(f"  No Sanger FASTA for {bc}, skipping.")
        continue
    if not os.path.isfile(bam_file):
        print(f"  No BAM for {bc}, skipping.")
        continue
    sanger_df = extract_sanger_codons(sanger_file[0])
    nano_df   = reconstruct_nanopore_codons_fast(bam_file)
    merged    = compare_and_save(bc, sanger_df, nano_df)
    all_merged.append(merged.assign(barcode=bc))
    plot_comparison(bc, merged)

# -----------------------------
# Summarize global comparisons
# -----------------------------
if not all_merged:
    print("No comparison data found!")
    exit()

global_df = pd.concat(all_merged, ignore_index=True)

# Save global
out_global = os.path.join(OUTPUT_DIR, "global_comparison_summary.csv")
global_df.to_csv(out_global, index=False)
print(f"Saved global summary -> {out_global}")

# Sanger-only variants
sanger_only = global_df[
    (global_df["Detected_by_Sanger"]) &
    (~global_df["Detected_by_Nanopore"])
].copy()

sanger_only["Difference(%)"] = sanger_only["frequency(%)_sanger"] - sanger_only["nano_freq(%)"]

out_sanger = os.path.join(OUTPUT_DIR, "sanger_only_variants.csv")
sanger_only.to_csv(out_sanger,index=False)
print(f"Sanger-only variants -> {out_sanger}")

# -----------------------------
# Combined Global Plot
# -----------------------------
plt.figure(figsize=(24,8))
global_df["label"] = global_df["barcode"] + ":" + \
                     global_df["position"].astype(str) + "_" + \
                     global_df["Sanger_codon"] + "/" + global_df["Nanopore_codon"]

plt.bar(global_df["label"], global_df["frequency(%)_sanger"], label="Sanger", alpha=0.7)
plt.bar(global_df["label"], global_df["nano_freq(%)"], label="Nanopore", alpha=0.5)

plt.xticks(rotation=90, fontsize=6)
plt.ylabel("Frequency (%)")
plt.title("Global Sanger vs Nanopore Codon Frequencies")
plt.legend()
plt.tight_layout()

out_global_png = os.path.join(OUTPUT_DIR, "global_codon_freq_comparison.png")
plt.savefig(out_global_png)
plt.close()
print(f"Saved global plot -> {out_global_png}")

# -----------------------------
# Summary Statistics
# -----------------------------
# Overall Pearson correlation
try:
    corr, pval = pearsonr(global_df["frequency(%)_sanger"],
                           global_df["nano_freq(%)"])
    print(f"Overall Pearson r: {corr:.3f} (p={pval:.2e})")
except Exception as e:
    print("Correlation failed:", e)

# Recovery rate
total_sanger = global_df["Detected_by_Sanger"].sum()
matched_sanger = ((global_df["Detected_by_Sanger"]) &
                  (global_df["Sanger_codon"] == global_df["Nanopore_codon"])).sum()
if total_sanger > 0:
    recovery_rate = matched_sanger/total_sanger*100
    print(f"Sanger → Nanopore recovery: {recovery_rate:.1f}% ({matched_sanger}/{total_sanger})")
else:
    print("No Sanger codons to compute recovery.")

print("Done.")
