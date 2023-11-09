from pyGeno.Genome import Genome  #Creating  personalized genome
from pyGeno.Gene import Gene #Printing all the proteins of a gene
from pyGeno.Transcript import Transcript
from pyGeno.Protein import Protein   #Printing all the proteins of a gene
from pyGeno.Exon import Exon

from pyGeno.SNPFiltering import SNPFilter  #filtering SNPs
from pyGeno.SNPFiltering import SequenceSNP  #Filtering SNPs
from pyGeno.tools.parsers.FastaTools import *

from argparse import ArgumentParser

def outDirCheck(parser, outDir) :
    outDir = os.path.abspath(outDir)
    if not os.path.isdir(outDir) :
        parser.error('The specified output directory %s does not exist!' %(outDir))
    else :
        return outDir

def get_parser():
    parser = ArgumentParser()
    parser.add_argument("-v",
                        dest="REF",
                        help="Version of the reference to be used (Default: GRCh38.78)",
                        default='GRCh38.78',
                        type=str)
    parser.add_argument("-s",
                        dest="SAMPLE",
                        help="Name of the sample (Default: Cancer)",
                        default='Cancer',
                        type=str)
    parser.add_argument("-snp",
                        dest="SNP_SET",
                        help="Name of the snp set of be used (Default: None)",
                        default=None)
    parser.add_argument("-a",
                        dest="ABUNDANCE",
                        help="relevant abundance.tsv file",
                        default='abundance.tsv',
                        type=str)
    parser.add_argument("-qual",
                        dest="SNP_QUALITY",
                        help="Minimun snp quality (Default: 20)",
                        default=20,
                        type=float)
    parser.add_argument("-tpm",
                        dest="TPM",
                        help="TPM threshold",
                        default=0,
                        type=float)
    parser.add_argument("-o",
                        dest="DIR_OUT",
                        help="Path to the output directory",
                        type=lambda x: outDirCheck(parser, x))
    return parser



class Qual_filter(SNPFilter):
    def __init__(self, threshold):
        super().__init__()
        self.threshold = threshold
    def filter(self, chromosome, **kwargs):
        for snp_set, snp in kwargs.items():
            if float(snp.quality) > self.threshold:
                return SequenceSNP(snp.alt.replace(',', ''))
        return None


def createPersoGenome(ref, snp, qual) :
    try :
        qual_filter = Qual_filter(qual)
        pG = Genome(name = ref, SNPs = snp, SNPFilter = qual_filter)
        #pG = Genome(name = ref, SNPs = snp)
        return pG
    except KeyError:
        return -1
    except ValueError:
        return -2



def main():

    global args
    args = get_parser().parse_args()
    print('Creating personalized genome...')
    persoGenome = createPersoGenome(args.REF, args.SNP_SET, args.SNP_QUALITY)
    if persoGenome == -1 :
        print( '    Genome %s is not installed in pyGeno' %(args.REF))
    elif persoGenome == -2 :
        print('    Snp set %s has not been imported in pyGeno' %(args.SNP_SET))
    else :
        print( '    Starting to export personalized proteome...')
        file = open(args.ABUNDANCE, "r")
        next(file)
        tpms = {}
        for line in file.readlines():
            a = line.strip().split("\t")
            enst = a[0].split(".")[0]
            exp = float(a[4])
            if exp > args.TPM :
                tpms[enst] = exp
        file.close()

        idk = 1
        outfile = args.DIR_OUT + '/' + "personalized_proteome.fasta"
        f = open(outfile,"a+")

        for trans in persoGenome.iterGet(Transcript):
            if trans.protein is None:
                continue
            elif trans.id in tpms:
                if not trans.exons:
                    continue
                exon_data = []
                for e in trans.exons:
                    if e is None:
                        continue
                    e_data = e.data
                    if not isinstance(e_data, list) or None in e_data:
                        continue
                    exon_data.append(e_data)
                uid = "".join([">",str(trans.chromosome.number),"_",str(idk)])
                fid = " ".join([uid, "Chromosome number:", trans.chromosome.number, "| Gene symbol:", trans.gene.name, ", Gene id:", trans.gene.id, "| Transcript id:", trans.id, "| Protein id:", trans.protein.id,"\n"])
                idk = idk + 1
                f.write(fid)
                f.write(trans.protein.sequence)
                f.write("\n")
        f.close()

main()
