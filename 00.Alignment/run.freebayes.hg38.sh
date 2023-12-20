sample=$1
bam=$2
outdir=$3

freebayes -f  homo_sapiens/hg38.fa -C 3 --min-coverage 10 $bam  >$outdir/$sample.vcf

perl 00.Alignment/convert.vcf.pl $sample $outdir/$sample.vcf $outdir/$sample.snps.txt $outdir/manifest.ini
cd $outdir
tar -zcvf $sample.snps.tar.gz $sample.snps.txt manifest.ini

python3 00.Alignment/import.snps.py $sample.snps.tar.gz $sample
