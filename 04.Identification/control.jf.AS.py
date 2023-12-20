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
import pysam
import random


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

def search_bam(protein_id, bam_file):
    matched_seqs = []
    bam = pysam.AlignmentFile(bam_file, "rb")
    bam_index = pysam.IndexedReads(bam)
    bam_index.build()
    for read in bam_index.find(protein_id):
            matched_seqs.append(read.query_sequence)
            print("matched seqs", read.query_sequence)
    return matched_seqs


def search_jf(query_seqs, jf_file):
    cmd = ['jellyfish', 'query', jf_file] + query_seqs
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

    output = result.stdout.strip()
    count = int(output.split()[1])
    if count > 0:
        return "yes"
    else:
        return "no"


def peptide_to_dna(peptide):
    dna_seqs = []
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]

    for aa in peptide:
        codons = [codon for codon, amino_acid in standard_table.forward_table.items() if amino_acid == aa]
        dna_seqs.append(codons)

    dna_seqs_combinations = list(itertools.product(*dna_seqs))
    dna_seqs_combinations = [''.join(dna_seq) for dna_seq in dna_seqs_combinations]

    for i, dna_seq_combination in enumerate(dna_seqs_combinations):
        if dna_seq_combination.startswith('ATG'):
            for j in range(len(dna_seq_combination) // 3):
                codon = dna_seq_combination[j*3:j*3+3]
                if codon in ['TAA', 'TAG', 'TGA']:
                    dna_seq_combination = dna_seq_combination[:j*3+3]
                    dna_seqs_combinations[i] = dna_seq_combination
                    break
        else:
            for j in range(len(dna_seq_combination) // 3):
                codon = dna_seq_combination[j*3:j*3+3]
                if codon == 'ATG':
                    dna_seq_combination = 'ATG' + dna_seq_combination[j*3:]
                    for k in range(len(dna_seq_combination) // 3):
                        codon = dna_seq_combination[k*3:k*3+3]
                        if codon in ['TAA', 'TAG', 'TGA']:
                            dna_seq_combination = dna_seq_combination[:k*3+3]
                            dna_seqs_combinations[i] = dna_seq_combination
                            break
                    break
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
    cancer_kmer_count = get_kmer_count(args.jf_cancer_file)

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
            protein_acs = [s for s in row['Protein AC'].split('/') if s]
            print("protein_acs", protein_acs)
            if len(protein_acs) > 20:
                protein_acs = protein_acs[:20]
                print("protein_acs_first20",protein_acs)

            for protein_ac in protein_acs:
                if protein_ac:
                    protein_id = protein_ac.rsplit("_", 1)[0]
                    print("protein_id", protein_id)

                    dna_seq = search_bam(protein_id, args.bam_file)
                    print("dna_seq",dna_seq)

                if dna_seq:
                    dna_seqs_combinations = peptide_to_dna(peptide)
                    dna_seq_change = dna_seq[:]
                    matched_dna = None
                    for dna_seq in dna_seq_change:
                        for dna_seq_combination in dna_seqs_combinations:
                            if dna_seq_combination in dna_seq:
                                matched_dna = dna_seq_combination
                                all_dna_seqs = ','.join(dna_seqs_combinations)
                                matched_control_jf = search_jf([matched_dna], args.jf_control_file)
                                matched_cancer_jf = search_jf([matched_dna], args.jf_cancer_file)
                                if matched_control_jf == 'yes':
                                    control_jf_count = calculate_count(matched_dna, args.jf_control_file)
                                    nFreq = count_norm(control_jf_count, control_kmer_count)
                                if matched_cancer_jf == 'yes':
                                    cancer_jf_count = calculate_count(matched_dna, args.jf_cancer_file)
                                    cFreq = count_norm(cancer_jf_count, cancer_kmer_count)

                                break

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
                            row['cFreq'] = cFreq
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
                        row['Matched in Cancer JF'] = 'NA'
                        row['TSA'] = 'NA'
                        row['nFreq'] = 'NA'
                        row['cFreq'] = 'NA'

                    writer.writerow(row)
                else:
                    print(f"Warning: Failed to find record with ID {protein_id} in the fasta file")

if __name__ == '__main__':
    main()
