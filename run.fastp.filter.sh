sample=$1
fq1=$2
fq2=$3
outdir=$4
#mkdir $outdir/$sample
cd $outdir
/mnt/disk2_workspace/jiangyiqi/bin/fastp -i $fq1 -o $outdir/$sample\_1.filter.fastq.gz -I $fq2 -O $outdir/$sample\_2.filter.fastq.gz -5 -3 -w 4
