#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name="dechat_all"
#SBATCH --mail-user=daniel.lule@lshtm.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6GB
#SBATCH --output=dechat_%j.out
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

module load anaconda
source activate dechat

INPUT_DIR="/home/lsh2202849/output_demultiplexed_nanoq"
OUTPUT_DIR="/home/lsh2202849/dechat_out"
THREADS=2
KMER=35

mkdir -p "$OUTPUT_DIR"

for INPUT_FASTQ in "$INPUT_DIR"/*.fastq "$INPUT_DIR"/*.fastq.gz; do
    [[ -e "$INPUT_FASTQ" ]] || continue

    base=$(basename "$INPUT_FASTQ")
    sample="${base%%.*}"

    SAMPLE_OUT_DIR="$OUTPUT_DIR/$sample"
    mkdir -p "$SAMPLE_OUT_DIR"

    echo "Processing $sample"

    cp "$INPUT_FASTQ" "$SAMPLE_OUT_DIR/input.fastq"

    cd "$SAMPLE_OUT_DIR"

    # run dechat INSIDE this directory
    dechat -i input.fastq -o "$sample" -t "$THREADS" -k "$KMER"

    cd -
done

echo "DONE"

