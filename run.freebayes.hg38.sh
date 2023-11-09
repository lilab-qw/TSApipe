sample=$1
bam=$2
outdir=$3

/home/wangmengyao/packages/freebayes/bin/freebayes -f  /mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/00.Alignment/homo_sapiens/hg38.fa -C 3 --min-coverage 10 $bam  >$outdir/$sample.vcf

perl /mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/00.Alignment/convert.vcf.pl $sample $outdir/$sample.vcf $outdir/$sample.snps.txt $outdir/manifest.ini
cd $outdir
tar -zcvf $sample.snps.tar.gz $sample.snps.txt manifest.ini

/home/wangmengyao/miniconda2/envs/test_py3/bin/python /mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/02.Canoical.cancer.proteome/import.snps.py $sample.snps.tar.gz $sample
