sample=$1
SRR=$2
outdir=/mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq
fq1=$outdir/$sample/${SRR}_1.fastq
fq2=$outdir/$sample/${SRR}_2.fastq
#outdir=/mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq
#mkdir $outdir/$sample

sh /mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/00.Alignment/run.fastp.filter.sh $sample $fq1 $fq2 $outdir/$sample && \
sh /mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/00.Alignment/run.star.hg38.sh $sample $outdir/$sample/$sample\_1.filter.fastq.gz $outdir/$sample/$sample\_2.filter.fastq.gz $outdir/$sample && \
sh /mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/00.Alignment/run.quant.sh $sample $outdir/$sample/$sample\_1.filter.fastq.gz $outdir/$sample/$sample\_2.filter.fastq.gz $outdir/$sample && \
sh /mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/00.Alignment/run.freebayes.hg38.sh $sample $outdir/$sample/$sample.bam $outdir/$sample && \
echo 'finsh $sample!'
