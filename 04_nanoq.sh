#!/bin/bash


# This Bash script trims Nanopore FASTQ files in a specified folder using the nanoq tool. 
# It processes both uncompressed (.fastq) and gzipped (.fastq.gz) files and outputs 
# trimmed FASTQs in the same directory. 

# Make the script executable, on the terminal run: chmod +x 04_nanoq.sh
# Run the script in the terminal: ./04_nanoq.sh


INPUT_DIR="/path/to/fastq" #input path to folder with fastq files

for f in "$INPUT_DIR"/*barcode*.fastq "$INPUT_DIR"/*barcode*.fastq.gz; do
    [ -e "$f" ] || continue

    base=$(basename "$f")
    base=${base%%.*}   # remove .fastq or .gz extension

    echo "Processing file: $f"

    nanoq \
        -i "$f" \
        -o "$INPUT_DIR/${base}.trimmed.fastq" \
        -q 15 \
        -l 500 \
        -m 1200
done

echo "All barcode FASTQ files processed."
