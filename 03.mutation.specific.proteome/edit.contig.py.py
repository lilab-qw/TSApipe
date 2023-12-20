
import sys
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(rna_editing_file, bam_file, output_fasta):
    rna_editing_sites = read_rna_editing_file(rna_editing_file)

    bam = pysam.AlignmentFile(bam_file, "rb")

    contigs = extract_contigs(rna_editing_sites, bam)

    SeqIO.write(contigs, output_fasta, "fasta")

def read_rna_editing_file(file_path):
    sites = []
    with open(file_path, "r") as f:
        next(f)  
        for line in f:
            fields = line.strip().split()
            site = (fields[0], int(fields[1]))
            sites.append(site)
    return sites

def extract_contigs(sites, bam, kmer_sizes=[90], min_mapq=0):
    contigs = []  

    for site in sites:
        chrom, pos = site

        for k in kmer_sizes:
            start, end = pos - k // 2, pos + k // 2
            reads = bam.fetch(chrom, start, end)

            for read in reads:
                if read.mapping_quality >= min_mapq:
                    seq = read.query_sequence
                    seq_length = len(seq)

                    contig_start = read.reference_start
                    contig_end = contig_start + seq_length
                    if contig_start <= start and contig_end >= end:
                        start_offset = start - contig_start
                        end_offset = end - contig_start
                        subseq = seq[start_offset:end_offset]
                        contig = SeqRecord(Seq(subseq), id=f"{chrom}:{start}-{end}_k{k}", description="")
                        print("contig",contig)
                        contigs.append(contig)

    return contigs

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_kmer_contigs.py <RNA_editing_file> <BAM_file> <output_fasta>")
        sys.exit(1)

    rna_editing_file = sys.argv[1]
    bam_file = sys.argv[2]
    output_fasta = sys.argv[3]

    main(rna_editing_file, bam_file, output_fasta)
