#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name="nanoHIV_demux"
#SBATCH --mail-user=daniel.lule@lshtm.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32GB
#SBATCH --output=dorado_demux_%j.out
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

module purge
module load dorado/0.9.1
module load samtools/1.20

echo "Dorado:   $(which dorado)"
echo "Samtools: $(which samtools)"

INPUT_DIR="/home/lsh2202849/hivdr"
OUTPUT_DIR="/home/lsh2202849/output_demultiplexed"
mkdir -p "$OUTPUT_DIR"

echo "Running Dorado demux…"
dorado demux \
    --kit-name SQK-NBD114-96 \
    --output-dir "$OUTPUT_DIR" \
    "$INPUT_DIR"

echo "Dorado demux completed."
echo "Converting BAM → FASTQ…"

# Handle both structures:
#   output/barcode01/*.bam
#   output/barcode01.bam
find "$OUTPUT_DIR" -type f -name "*.bam" | while read -r bam; do
    echo "Found BAM: $bam"

    barcode=$(basename "$bam" .bam)
    folder=$(dirname "$bam")
    out_fastq="$folder/${barcode}.fastq"

    echo "Converting → $out_fastq"

    # -f 4 ensures unmapped reads are output (dorado demux BAM files are unmapped)
    samtools fastq -f 4 -n "$bam" > "$out_fastq"

    if [[ ! -s "$out_fastq" ]]; then
        echo "⚠ WARNING: FASTQ is empty for $barcode — BAM likely contains no reads."
    fi
done

echo "All BAM→FASTQ conversions completed."

