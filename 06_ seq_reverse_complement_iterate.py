#!/usr/bin/env python3

"""
This script processes one or more FASTA files in a given directory and checks each 
sequence against a reference sequence. Sequences that align better in the 
reverse-complement orientation are automatically reversed and complemented. The 
oriented sequences are written to new FASTA files in an output directory. 

Example usage:

python seq_reverse_complement_multi.py \
  -i input_fastas/ \
  -r hxb2_reference.fasta \
  -o oriented_fastas/
"""

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


# --- Helper function ---
def is_reverse(seq, reference, aligner):
    original_score = aligner.score(seq, reference)
    revcomp_seq = str(Seq(seq).reverse_complement())
    revcomp_score = aligner.score(revcomp_seq, reference)
    return revcomp_score > original_score


def main():

    # --- Parse arguments ---
    parser = argparse.ArgumentParser(
        description="Orient sequences in multiple FASTA files to a single reference."
    )
    parser.add_argument(
        "-i", "--input_dir", required=True,
        help="Directory containing input FASTA files"
    )
    parser.add_argument(
        "-r", "--reference", required=True,
        help="Reference FASTA file (e.g. HXB2)"
    )
    parser.add_argument(
        "-o", "--output_dir", required=True,
        help="Directory for oriented output FASTA files"
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # --- Load reference once ---
    ref_record = SeqIO.read(args.reference, "fasta")
    ref_seq = str(ref_record.seq)

    # --- Initialize aligner once ---
    aligner = PairwiseAligner()
    aligner.mode = "global"

    # --- Iterate over FASTA files ---
    for fasta_file in sorted(os.listdir(args.input_dir)):
        if not fasta_file.lower().endswith((".fa", ".fasta", ".fna")):
            continue

        input_path = os.path.join(args.input_dir, fasta_file)
        output_path = os.path.join(
            args.output_dir,
            fasta_file.replace(".fasta", "_oriented.fasta")
        )

        print(f"\nProcessing {fasta_file}")

        corrected_records = []

        for record in SeqIO.parse(input_path, "fasta"):
            seq = str(record.seq)

            if is_reverse(seq, ref_seq, aligner):
                print(f"  {record.id} is reverse → reversing")
                record.seq = record.seq.reverse_complement()

            corrected_records.append(record)

        # --- Write output (single-line FASTA) ---
        with open(output_path, "w") as out:
            for record in corrected_records:
                out.write(f">{record.id}\n{record.seq}\n")

        print(f"  Saved → {output_path}")

    print("\nAll FASTA files processed successfully.")


if __name__ == "__main__":
    main()
