#!/bin/bash

# This Bash script aligns per-barcode read FASTA files to their corresponding reference 
# FASTA sequences using minimap2, producing sorted and indexed BAM files. The script 
# automatically handles reference indexing and processes all matching barcode references 
# in a specified directory.

set -euo pipefail
shopt -s nullglob nocaseglob  # ensures pattern expansion and case-insensitive matching

REF_DIR="fasta_pr_references_medaka" # set path to Sanger references directory oriented as forward reads
READ_DIR="output_dechat" # set path to input directory of Nanopore reads
OUT_DIR="minimap2_results2" # set path to output directory

THREADS=8
mkdir -p "$OUT_DIR"

# Loop over all reference fasta files matching the pattern, case-insensitive
for REF in "$REF_DIR"/barcode*_ref_*.fasta; do
    # Skip if no files matched
    [[ -e "$REF" ]] || continue

    # Extract barcode prefix before the first '_ref_'
    BASENAME=$(basename "$REF")
    BARCODE="${BASENAME%%_ref_*}"

    READ_FASTA="${READ_DIR}/${BARCODE}/${BARCODE}.fasta"

    echo "Processing $BARCODE"

    if [[ ! -f "$READ_FASTA" ]]; then
        echo "  ❌ Reads not found: $READ_FASTA"
        continue
    fi

    # -------------------------
    # Index reference if needed
    # -------------------------
    if [[ ! -f "${REF}.fai" ]]; then
        echo "  Indexing reference FASTA"
        samtools faidx "$REF"
    fi

    # Optional minimap2 index (speed-up)
    if [[ ! -f "${REF}.mmi" ]]; then
        minimap2 -d "${REF}.mmi" "$REF"
    fi

    # -------------------------
    # Alignment → sorted BAM
    # -------------------------
    BAM_OUT="${OUT_DIR}/${BARCODE}.sorted.bam"

    minimap2 -ax map-ont "${REF}.mmi" "$READ_FASTA" \
        | samtools sort -@ "$THREADS" -o "$BAM_OUT"

    samtools index "$BAM_OUT"

    echo "  ✅ Output: $BAM_OUT"
    echo
done

