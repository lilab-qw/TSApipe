import sys
import csv
import re
import os
from Bio import SeqIO
import argparse

def main(csv_path, fasta_path, output_path):
    csv.field_size_limit(sys.maxsize)
    for path in [csv_path, fasta_path]:
        if not os.path.isfile(path):
            print(f"Error: '{path}' file not found.")
            sys.exit(1)
    sequences = []
    with open(fasta_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(str(record.seq))
    sequences = list(set(sequences))

    with open(csv_path, 'r') as csv_file, open(output_path, 'w') as output_file:
        reader = csv.DictReader(csv_file, delimiter='\t')
        fieldnames = reader.fieldnames

        if "Sq" not in fieldnames:
            print("Error: 'Sq' column not found in CSV file.")
            sys.exit(1)

        writer = csv.DictWriter(output_file, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in reader:
            peptide = row["Sq"]
            is_peptide_found = False
            for seq in sequences:
                if peptide in seq:
                    is_peptide_found = True
                    break
            if not is_peptide_found:
                print(f"{peptide} is not in any sequence")
                writer.writerow(row)
            else:
                print(f"{peptide} is in at least one sequence: {seq}")
                    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter CSV rows by matching peptides in a FASTA file.')
    parser.add_argument('csv_path', metavar='<CSV file>', help='Path to CSV file')
    parser.add_argument('fasta_path', metavar='<FASTA file>', help='Path to FASTA file')
    parser.add_argument('output_path', metavar='<output file>', help='Path to output file')
    args = parser.parse_args()

    main(args.csv_path, args.fasta_path, args.output_path)
