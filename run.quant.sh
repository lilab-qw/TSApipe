sample=$1
fq1=$2
fq2=$3
outdir=$4
/mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/software/kallisto/kallisto  quant -i /mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/00.Alignment/homo_sapiens/transcriptome.idx -o $outdir/ -t 4 $fq1 $fq2
