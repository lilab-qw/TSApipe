sample=$1
outdir=$3
fq1=$2
fq2=$3
mkdir $outdir/$sample

sh 00.Alignment/run.fastp.filter.sh $sample $fq1 $fq2 $outdir/$sample && \
sh 00.Alignment/run.star.hg38.sh $sample $outdir/$sample/$sample\_1.filter.fastq.gz $outdir/$sample/$sample\_2.filter.fastq.gz $outdir/$sample && \
sh 00.Alignment/run.quant.sh $sample $outdir/$sample/$sample\_1.filter.fastq.gz $outdir/$sample/$sample\_2.filter.fastq.gz $outdir/$sample && \
echo 'finsh $sample!'
