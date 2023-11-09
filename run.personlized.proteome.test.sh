sample=$1
snp=$2
outdir=/mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/$sample/Canoical.cancer.proteome
abundance=/mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/$sample/abundance.tsv
cd /mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/$sample
python /mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/02.Canoical.cancer.proteome/personlized.proteome.test.py -v GRCh38.78   -s $sample -o $outdir -snp $snp -a $abundance  -tpm 0 -qual 20
echo $(date) finished..
