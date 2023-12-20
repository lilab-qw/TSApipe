sample=$1
bam=$2

freebayes -f  homo_sapiens/hg38.fa -C 3 --min-coverage 10 $bam  > $outdir/$sample.vcf

python snp.contig.py  $outdir/$sample.vcf $bam  ${sample}.snp.contig

echo $(date) ${sample}.contig finished ...

awk '{if (/^>/) {if (seq) {print name"\n"seq}; name=$0; seq=""; i++} else {seq = seq $0}} END {if (seq) {print name"\n"seq}}' ${sample}.snp.contig | awk '!a[$0]++' | awk '{if (/^>/) {name=$0} else {print name"\n"$0}}' > ${sample}.snp.contig.fasta

echo $(date) ${sample}.contig duplication finished ...

EMBOSS-6.6.0/emboss/transeq  ${sample}.snp.contig.fasta  -frame=3 ${sample}.snp.fasta