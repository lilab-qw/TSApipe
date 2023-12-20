sample=$1
outdir=$2
AS_type=$3
LCOUNT=$4

python 04.Identification/control.fasta.py  ${outdir}/${sample}.capaMHC.MHCbinding.final.csv Control/annotationFiles/control_control_reformat_70char.fasta ${outdir}/${sample}.control.fasta.csv

python 04.Identification/control.jf.edit-snp-indel.py ${outdir}/rna-edit/${sample}.control.fasta.csv Control/control.trim.24.jf ${outdir}/${sample}.trim.24.jf   ${outdir}/rna-edit/${sample}.control.jf.csv ${outdir}/rna-edit/${sample}.candidate.contig.fasta

python 04.Identification/control.jf.edit-snp-indel.py $outdir/indel/${sample}.control.fasta.csv Control/control.trim.24.jf $outdir/${sample}.trim.24.jf   $outdir/indel/${sample}.control.jf.csv $outdir/indel/${sample}.indel.contig.fasta

python 04.Identification/control.jf.edit-snp-indel.py $outdir/snp/${sample}.control.fasta.csv Control/control.trim.24.jf $outdir/${sample}.trim.24.jf   $outdir/snp/${sample}.control.jf.csv $outdir/snp/${sample}.snp.contig.fasta

python 04.Identification/control.jf.fusion.py $outdir/GENE-FUSION/${sample}.control.fasta.csv Control/control.trim.24.jf $outdir/${sample}.trim.24.jf   $outdir/GENE-FUSION/${sample}.control.jf.csv $outdir/star-fusion/FusionInspector/${sample}.finspector.fa

python 04.Identification/control.jf.AS.py $outdir/rmats-${AS_type}/${sample}.control.fasta.csv Control/control.trim.24.jf $outdir/${sample}.trim.24.jf   $outdir/rmats-${AS_type}/${sample}.control.jf-3.csv $outdir/${sample}.bam

python 04.Identification/control.jf.spe.py $outdir/spe/${sample}.control.fasta.csv Control/control.trim.24.jf $outdir/${sample}.trim.24.jf   $outdir/spe/${sample}.control.jf.csv $outdir/linear_${sample}.${LCOUNT}.control.0/assembly.fasta
