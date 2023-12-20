from Bio.Seq import Seq
import argparse
import csv
import re
import subprocess
import sys
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio import SeqIO
import tempfile
import os
import itertools

csv.field_size_limit(sys.maxsize) 
arser = argparse.ArgumentParser(description='Search for peptides in jellyfish class files')
parser.add_argument('csv_file', help='Input CSV file')
parser.add_argument('jf_control_file', help='Input control jellyfish class file')
parser.add_argument('jf_cancer_file', help='Input cancer jellyfish class file')
parser.add_argument('result_file', help='Output CSV file for results')
parser.add_argument('fasta_file', help='Input fasta file')
args = parser.parse_args()

def parse_csv(csv_file):
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f,delimiter='\t')
        data = list(reader)
    return data


def search_fasta(protein_id, fasta_file):
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id.startswith(protein_id):
            return str(record.seq)
    return None

def search_jf(query_seqs, jf_file):
    cmd = ['jellyfish', 'query', jf_file] + query_seqs
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    output = result.stdout.strip()
    count = int(output.split()[1])
    print('search_jf_count',count)
    if count > 0:
        return "yes"
    else:
        return "no"

def peptide_to_dna(peptide):
    dna_seqs = []
    standard_table = CodonTable.unambiguous_dna_by_id[1]

    for aa in peptide:
        codons = [codon for codon, amino_acid in standard_table.forward_table.items() if amino_acid == aa]
        dna_seqs.append(codons)

    dna_seqs_combinations = list(itertools.product(*dna_seqs))
    dna_seqs_combinations = [''.join(dna_seq) for dna_seq in dna_seqs_combinations]
    return dna_seqs_combinations


def get_kmer_count(jf_file):
    cmd = ['jellyfish', 'stats', jf_file]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    output = result.stdout.strip()
    kmer_count = int(output.split()[1])

    return kmer_count

def calculate_count(dna_seq, jf_file):
    cmd = ['jellyfish', 'query', jf_file, dna_seq]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    output = result.stdout.strip()
    print('output',output)
    count = int(output.split()[1])

    return count

def count_norm(count, total):
    res = (float(count) / total) * 10**9
    return res

def main():
    if not os.path.isfile(args.csv_file):
        print(f"Error: Input CSV file does not exist or is not a file: {args.csv_file}")
        sys.exit(1)
    if not os.path.isfile(args.jf_control_file):
        print(f"Error: Input control jellyfish class file does not exist or is not a file: {args.jf_control_file}")
        sys.exit(1)
    if not os.path.isfile(args.jf_cancer_file):
        print(f"Error: Input cancer jellyfish class file does not exist or is not a file: {args.jf_cancer_file}")
        sys.exit(1)
    if not os.path.isfile(args.fasta_file):
        print(f"Error: Input fasta file does not exist or is not a file: {args.fasta_file}")
        sys.exit(1)

    control_kmer_count = get_kmer_count(args.jf_control_file)
    print('control_kmer_count',control_kmer_count)
    cancer_kmer_count = get_kmer_count(args.jf_cancer_file)
    print("cancer_kmer_count",cancer_kmer_count)
    csv_data = parse_csv(args.csv_file)

    with open(args.result_file, 'w', newline='') as f:
        fieldnames = csv_data[0].keys()
        fieldnames = list(fieldnames) + ['Matched DNA', 'Matched in Control JF', 'Matched in Cancer JF', 'nFreq', 'cFreq', 'TSA']
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for row in csv_data:
            peptide = row['Sq']
            print("peptide",peptide)
            protein_id = row['Protein AC']
            protein_id = protein_id.split("_")[0]
            print("protein_id",protein_id)

            dna_seq = search_fasta(protein_id, args.fasta_file)
            print("dna_seq",dna_seq)

            if dna_seq:
                dna_seqs_combinations = peptide_to_dna(peptide)
                print("dna_seqs_combinations",dna_seqs_combinations)

                matched_dna = None
                for dna_seq_combination in dna_seqs_combinations:
                    if dna_seq_combination in dna_seq:
                        matched_dna = dna_seq_combination
                        all_dna_seqs = ','.join(dna_seqs_combinations)
                        matched_control_jf = search_jf([matched_dna], args.jf_control_file)
                        print('matched_control_jf',matched_control_jf)
                        matched_cancer_jf = search_jf([matched_dna], args.jf_cancer_file)
                        print('matched_cancer_jf',matched_cancer_jf)
                        if matched_control_jf == 'yes':
                            control_jf_count = calculate_count(matched_dna, args.jf_control_file)
                            print('control_jf_count',control_jf_count)
                            nFreq = count_norm(control_jf_count, control_kmer_count)
                        if matched_cancer_jf == 'yes':
                            cancer_jf_count = calculate_count(matched_dna, args.jf_cancer_file)
                            print('cancer_jf_count',cancer_jf_count)
                            cFreq = count_norm(cancer_jf_count, cancer_kmer_count)
                        break
                    else:
                        print("Warning: Failed to find a possible sequence in the FASTA file")

                if matched_dna:

                    if matched_control_jf == 'yes':
                        if matched_cancer_jf == 'yes':
                            row['Matched in Control JF'] = 'y'
                            row['Matched in Cancer JF'] = 'y'
                            row['TSA'] = 'maybeyes'
                            row['nFreq'] = nFreq
                            row['cFreq'] = cFreq
                        else:
                            row['Matched in Control JF'] = 'y'
                            row['Matched in Cancer JF'] = 'n'
                            row['TSA'] = 'NA'
                            row['nFreq'] = nFreq
                            row['cFreq'] = 'NA'
                    elif matched_cancer_jf == 'yes':
                        row['Matched in Control JF'] = 'n'
                        row['Matched in Cancer JF'] = 'y'
                        row['TSA'] = 'yes'
                        row['cFreq'] = nFreq
                        row['nFreq'] = 'NA'
                    else:
                        row['Matched in Control JF'] = 'n'
                        row['Matched in Cancer JF'] = 'n'
                        row['TSA'] = 'NA'
                        row['nFreq'] = 'NA'
                        row['cFreq'] = 'NA'
                    row['Matched DNA'] = matched_dna

                else:
                    row['Matched DNA'] = 'NA'
                    row['Matched in Control JF'] = 'NA'
                    row['nFreq'] = 'NA'
                    row['cFreq'] = 'NA'
                    row['Matched in Cancer JF'] = 'NA'
                    row['TSA'] = 'NA'

                writer.writerow(row)
            else:
                print(f"Warning: Failed to find record with ID {protein_id} in the fasta file")

if __name__ == '__main__':
    main()
