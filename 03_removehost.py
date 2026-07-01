#!/usr/bin/env python

'''
This script filters host reads from Nanopore sequencing FASTQ files by 
aligning them to a reference genome (human genome "human_g1k_v37.fasta") and 
extracting unmapped reads. It processes all FASTQ files in a specified input folder.


Usage python 03_removehost.py -i /path/to/input_folder -r human_g1k_v37.fasta -o /path/to/output_folder

'''

import sys
import os
import glob
import argparse

def main():
    cli = argparse.ArgumentParser()

    cli.add_argument('-i', '--InputFolder', help="Folder containing basecalled Nanopore fastq files. Only files ending in .fastq will be used", required=True)
    cli.add_argument('-r', '--Reference', help="Host Reference fasta or fasta.gz file", required=True)
    cli.add_argument('-o', '--OutputFolder', help="Output Folder. Default is ~/dehost_output/test", required=False, default='~/dehost_output/test')

    cli.add_argument('-t', '--threads', help="Number of threads. Default is 4", type=int, required=False, default=4)
    method = cli.add_mutually_exclusive_group()
    method.add_argument('--Nanopore', help="Select if you used Nanopore Sequencing", action='store_const', dest='seq_method', const='map-ont', default='map-ont')
    method.add_argument('--PacBio', help="Select if you used PacBio Sequencing", action='store_const', dest='seq_method', const='map-pb')
    args = cli.parse_args()

    # Expand user path
    output_folder = os.path.expanduser(args.OutputFolder)
    os.makedirs(output_folder, exist_ok=True)

    files = glob.glob(os.path.join(args.InputFolder, "*.fastq"))

    for fastq_file in files:
        base = os.path.splitext(os.path.basename(fastq_file))[0]
        sam_path = os.path.join(output_folder, f"{base}.sam")
        unmapped_ids = os.path.join(output_folder, f"{base}_unmapped_ids.txt")
        final_fastq = os.path.join(output_folder, f"{base}_filtered.fastq")

        # Step 1: Align with minimap2
        minimap2_cmd = f"minimap2 -ax {args.seq_method} -t {args.threads} {args.Reference} {fastq_file} > {sam_path}"
        print(f"[INFO] Running minimap2: {minimap2_cmd}")
        os.system(minimap2_cmd)

        # Step 2: Extract unmapped read IDs
        extract_ids_cmd = f"samtools view -f 4 {sam_path} | cut -f1 | sort | uniq > {unmapped_ids}"
        print(f"[INFO] Extracting unmapped read IDs: {extract_ids_cmd}")
        os.system(extract_ids_cmd)

        # Step 3: Extract full reads from original FASTQ using seqtk
        extract_fastq_cmd = f"seqtk subseq {fastq_file} {unmapped_ids} > {final_fastq}"
        print(f"[INFO] Extracting full reads with original headers: {extract_fastq_cmd}")
        os.system(extract_fastq_cmd)

        # Optional: delete intermediate SAM file
        os.remove(sam_path)
        os.remove(unmapped_ids)

        print(f"[DONE] Filtered FASTQ written to: {final_fastq}")

if __name__ == "__main__":
    sys.exit(main())
