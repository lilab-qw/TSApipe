sample=$1
outdir=$2
bam=$3
python reditools2.0/src/cineca/reditools.py -f  $bam  -r Homo_sapiens/hg38/hg38.fa  -o ${outdir}/rna-edit/${sample}.rnaedit

REDItools/bin/selectPositions.py -i ${outdir}/rna-edit/${sample}.rnaedit -r -c 10 -v 1 -f 0.1 -o ${outdir}/rna-edit/${sample}.candidate.txt

python edit.contig.py ${outdir}/rna-edit/${sample}.candidate.txt $bam  ${outdir}/rna-edit/${sample}.candidate.contig
awk '{if (/^>/) {if (seq) {print name"\n"seq}; name=$0; seq=""; i++} else {seq = seq $0}} END {if (seq) {print name"\n"seq}}'  ${outdir}/rna-edit/${sample}.candidate.contig | awk '!a[$0]++' | awk '{if (/^>/) {name=$0} else {print name"\n"$0}}' >  ${outdir}/rna-edit/${sample}.candidate.contig.fasta
EMBOSS-6.6.0/emboss/transeq  ${outdir}/rna-edit/${sample}.candidate.contig.fasta  -frame=3 ${outdir}/rna-edit/${sample}.candidate.fasta
echo $(date) ${sample}_candidadate_fasta finished

