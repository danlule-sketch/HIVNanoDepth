#!/bin/bash
set -euo pipefail

#-----------------------------------------------------------------------------------------
# Systematic BAM subsampling
#-----------------------------------------------------------------------------------------
# This Bash script performs systematic subsampling of Nanopore BAM files
# across multiple predefined sampling proportions.


# Input BAM files are read from: minimap2_results2_medaka/

# The script:

# 1. Iterates through all Barcode*.sorted.bam files
# 2. Subsamples reads at multiple sequencing-depth proportions using samtools view
# 3. Generates reproducible subsamples using a fixed random seed
# 4. Indexes each subsampled BAM file using samtools index

# The resulting BAM files can be used for:

# sequencing depth sensitivity analyses
# benchmarking variant-calling robustness
# evaluating performance stability under reduced read depth

#-----------------------------------------------------------------------------------------


BAM_DIR=/Users/daniellulebugembe/Documents/PhD/reading7_simulation/minimap2_results2_medaka

FRACTIONS="0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.4 0.5"

OUTDIR=subsampled_bams
mkdir -p "$OUTDIR"

for bam in "$BAM_DIR"/Barcode*.sorted.bam; do

    bc=$(basename "$bam" .sorted.bam)

    echo "Processing $bc"

    for frac in $FRACTIONS; do

        safe_frac=$(echo "$frac" | tr '.' 'p')

        outbam="${OUTDIR}/${bc}_${safe_frac}.bam"

        echo "  Subsampling at $frac"

        samtools view \
            -@ 4 \
            -bs 42${frac#0} \
            "$bam" \
            -o "$outbam"

        samtools index "$outbam"

    done

done

echo "Done."