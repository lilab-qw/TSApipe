import sys

from pyGeno.importation.SNPs import *
#snpdatawrap="/mnt/disk2_workspace/wangmengyao/wgTSA/test/Tumor/Tumor.snps.tar.gz"
#snpdatawrap="/mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/mammary_result/Tumor/Tumor.snps.tar.gz"
snpdatawrap = sys.argv[1]
#snps = sys.argv[2]
#deleteSNPs(snps)
importSNPs(snpdatawrap)

# the database wiil be locked while import snps, so no parall 

