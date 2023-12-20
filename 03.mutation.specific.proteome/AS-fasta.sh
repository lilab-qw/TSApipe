sample=$1
AS_type=$2
bam=$3
outdir=$4
mkdir rmats
/disk2/apps/software/Python/2.7.18-GCCcore-9.3.0/bin/python  rMATS.4.0.1/rMATS-turbo-Linux-UCS4/rmats.py --b1 ${sample}.txt --b2 control.txt --gtf homo_sapiens/Homo_sapiens.GRCh38.96.gtf  --od ${outdir}/rmats  -t paired --nthread 7   --anchorLength 1 --tstat 6
awk -F '\t' '$20 == "FDR" || $20<=0.05 {print}' ${outdir}/rmats/${AS_type}.MATS.JC.txt > $outdir/rmats/${sample}.${AS_type}.MATS.JC.filter.txt
cut -f 4,6,7 ${outdir}/rmats/${sample}.${AS_type}.MATS.JC.filter.txt > ${AS_type}.out.bed
sed -i '1d' ${AS_type}.out.bed
bedtools intersect -a $bam  -b ${AS_type}.out.bed > ${AS_type}.bam
samtools view ${AS_type}.bam -O SAM |  awk -F"\t" '{print ">"$1"\n"$10}' > ${AS_type}.bam.target.fasta
EMBOSS-6.6.0/emboss/transeq ${AS_type}.bam.target.fasta -frame=3 ${sample}.${AS_type}.fasta
rm -r ${AS_type}.out.bed  ${AS_type}.bam  ${AS_type}.bam.target.fasta
