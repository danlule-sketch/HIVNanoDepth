#!/usr/bin/env bash

# This Bash script performs per-barcode variant calling on Nanopore sequencing data 
# using Clair3. It processes sorted BAM files, matches each to its corresponding reference 
# FASTA, and generates variant calls using a pre-trained Clair3 model optimized for 
# ONT R10.4.1 SUP data. The Clair3 model used is in a folder named "/clair3_models"

# ------------------------------
# Activate Clair3 conda environment
# ------------------------------
CONDA_ENV="clair3"
source ~/anaconda3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"

set -euo pipefail
shopt -s nocaseglob nullglob

# ------------------------------
# Directories (HPC paths)
# ------------------------------
BAM_DIR="/home/lsh2202849/minimap2_results2"  # input Minimap2 output BAM file directory
REF_DIR="/home/lsh2202849/fasta_pr_references_medaka" # input Medaka consensus reference 
OUT_DIR="/home/lsh2202849/clair3_variants"
MODEL_DIR="$HOME/clair3_models/r1041_e82_400bps_sup_v500"

THREADS=8
mkdir -p "$OUT_DIR"

# ------------------------------
# Main loop
# ------------------------------
for BAM in "${BAM_DIR}"/barcode*.sorted.bam; do
    [[ -e "$BAM" ]] || continue

    # Extract barcode (barcode01)
    BARCODE=$(basename "$BAM" .sorted.bam)

    # Capitalize first letter  Barcode01
    BARCODE_CAP="$(tr '[:lower:]' '[:upper:]' <<< "${BARCODE:0:1}")${BARCODE:1}"

    # Find matching FASTA only (ignore .fai/.mmi)
    REF_MATCH=("${REF_DIR}/${BARCODE_CAP}_ref_"*.fasta)

    if [[ ${#REF_MATCH[@]} -eq 0 ]]; then
        echo " No reference found for $BARCODE"
        continue
    fi

    REF="${REF_MATCH[0]}"
    OUT="${OUT_DIR}/${BARCODE}"

    echo "===================================="
    echo "Running Clair3 for $BARCODE"
    echo "  BAM: $BAM"
    echo "  REF: $REF"
    echo "  OUT: $OUT"
    echo "===================================="

    run_clair3.sh \
    --bam_fn "$BAM" \
    --ref_fn "$REF" \
    --threads "$THREADS" \
    --platform ont \
    --model_path "$MODEL_DIR" \
    --output "$OUT" \
    --haploid_precise \
    --include_all_ctgs

done

