sample=$1
fq1=$2
fq2=$3
dir=$4
outdir=$dir/$sample
mkdir $outdir

/disk2/apps/software/Python/2.7.18-GCCcore-9.3.0/bin/python  /home/wangmengyao/miniconda2/bin/OptiTypePipeline.py -i $fq1 $fq2 --rna --outdir $outdir --prefix $sample
