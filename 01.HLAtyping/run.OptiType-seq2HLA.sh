sample=$1
fq1=$2
fq2=$3
outdir=$4
mkdir $outdir

python2  OptiTypePipeline.py -i $fq1 $fq2 --rna --outdir $outdir --prefix $sample
python3  seq2HLA.py  -1 $fq1 -2 $fq2  -r ${sample} -p 3
