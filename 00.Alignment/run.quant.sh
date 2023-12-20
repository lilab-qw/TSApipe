sample=$1
fq1=$2
fq2=$3
outdir=$4
kallisto  quant -i /homo_sapiens/transcriptome.idx -o $outdir/ -t 4 $fq1 $fq2
