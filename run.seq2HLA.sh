sample=$1
fq1=/mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/${sample}/${sample}_1.filter.fastq.gz
fq2=/mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/${sample}/${sample}_2.filter.fastq.gz
python /mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/software/seq2HLA/seq2HLA.py  -1 $fq1 -2 $fq2   -r ${sample} -p 3
