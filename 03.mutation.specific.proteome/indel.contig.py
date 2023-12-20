import sys
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(vcf_file, bam_file, output_fasta):
    rna_editing_sites = read_vcf_file(vcf_file)

    bam = pysam.AlignmentFile(bam_file, "rb")

    contigs = extract_contigs(rna_editing_sites, bam)

    SeqIO.write(contigs, output_fasta, "fasta")

def read_vcf_file(file_path):
    sites = []
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue  
            fields = line.strip().split()
            info = parse_info(fields[7])
            print("info",info)
            if info["TYPE"] == "ins" or  info["TYPE"] == "del" or info["TYPE"] == "complex" :
                print("info[type]",info["TYPE"])
                site = (fields[0], int(fields[1]))
                print("site",site)
                sites.append(site)
    return sites

def parse_info(info_str):
    info = {}
    for item in info_str.split(";"):
        key, value = item.split("=")
        info[key] = value
    return info

def extract_contigs(sites, bam, kmer_sizes=[90], min_mapq=20):
    contigs = []  
    for site in sites:
        chrom, pos = site

        for k in kmer_sizes:
            start, end = pos - k // 2, pos + k // 2

            for pileupcolumn in bam.pileup(chrom, start, end, max_depth=1000000):
                for pileupread in pileupcolumn.pileups:
                    if pileupread.query_position is not None:
                        read = pileupread.alignment
                        if read.mapping_quality >= min_mapq:
                            seq = read.query_sequence
                            seq_length = len(seq)

                            contig_id = f"{read.query_name}_k{k}_{read.is_reverse}"
                            
                            contig_start = read.reference_start
                            contig_end = contig_start + seq_length
                            if contig_start <= start and contig_end >= end:
                               
                                start_offset = start - contig_start
                                end_offset = end - contig_start
                                subseq = seq[start_offset:end_offset]
                                contig = SeqRecord(Seq(subseq), id=f"{chrom}:{start}-{end}_k{k}", description="")
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
