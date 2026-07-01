#!/bin/bash
set -euo pipefail

# This Bash script generates per-barcode consensus sequences from Nanopore sequencing 
# data using Medaka. It processes sorted BAM files, matches each to its corresponding 
# reference FASTA, and runs Medaka consensus in a controlled Conda environment.

# -------------------------------
# Activate Conda environment
# -------------------------------
source /home/lsh2202849/anaconda3/etc/profile.d/conda.sh
conda activate medaka

# -------------------------------
# Paths and parameters
# -------------------------------
BAM_DIR="/home/lsh2202849/bam" # input your BAM file directory
REF_DIR="/home/lsh2202849/medaka_ref" # input your Sanger sequence reference directory  
OUT_DIR="/home/lsh2202849/medaka_results" # input your output directory
THREADS=2  # Reduced threads for CPU-only environment
MODEL="r104_e81_sup_g5015"

mkdir -p "$OUT_DIR"

# -------------------------------
# Get BAM list
# -------------------------------
mapfile -t BAM_LIST < <(ls -1 ${BAM_DIR}/barcode*.sorted.bam ${BAM_DIR}/barcode*.sorted_by_coord.bam | sort)
NUM_BAMS=${#BAM_LIST[@]}

# -------------------------------
# Debugging BAM List and Barcode Extraction
# -------------------------------
echo "Debug: BAM files found:"
for BAM in "${BAM_LIST[@]}"; do
    echo "Full path BAM file: $BAM"
    
    # Extract barcode by removing ".sorted.bam" and ".sorted_by_coord.bam" suffixes
    BARCODE=$(basename "$BAM" | sed -e 's/.sorted_by_coord.bam//' -e 's/.sorted.bam//')
    echo "Extracted Barcode: $BARCODE"

    if [[ -z "$BARCODE" ]]; then
        echo "ERROR: Barcode extraction failed for $BAM"
        continue
    fi

    # -------------------------------
    # Match reference
    # -------------------------------
    REF=$(ls ${REF_DIR}/Barcode${BARCODE#barcode}_ref_*.fasta 2>/dev/null)
    echo "Looking for reference in $REF_DIR for barcode $BARCODE"

    if [[ -n "$REF" ]]; then
        echo "Reference found: $REF"
    else
        echo "ERROR: Reference not found for $BARCODE"
        continue
    fi

    # -------------------------------
    # Setup output directory
    # -------------------------------
    OUT="${OUT_DIR}/${BARCODE}"
    mkdir -p "$OUT"
    echo "Output directory created at: $OUT"

    # -------------------------------
    # Coordinate-sort BAM (if needed)
    # -------------------------------
    COORD_BAM="${BAM%.sorted_by_coord.bam}.sorted_by_coord.bam"
    if [[ ! -f "$COORD_BAM" ]]; then
        samtools sort -@ "$THREADS" -o "$COORD_BAM" "$BAM"
        samtools index "$COORD_BAM"
    fi

    # -------------------------------
    # Run Medaka consensus
    # -------------------------------
    echo "Running Medaka consensus"
    echo "  Barcode : $BARCODE"
    echo "  BAM     : $COORD_BAM"
    echo "  REF     : $REF"
    echo "  OUT     : $OUT"
    echo "  THREADS : $THREADS"

    MEDAKA_OUT="${OUT}/${BARCODE}_consensus.h5"
    medaka consensus "$COORD_BAM" "$MEDAKA_OUT" --threads $THREADS --model $MODEL

    echo "Consensus file saved as: $MEDAKA_OUT"

done

