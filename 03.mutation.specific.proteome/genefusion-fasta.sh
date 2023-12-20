sample=$1
outdir=$2
fq1=$3
fq2=$4

STAR-Fusion-v1.10.0/STAR-Fusion   --genome_lib_dir homo_sapiens/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir \
             --left_fq  $fq1 \
             --right_fq $fq2 \
             --output_dir  $outdir/star-fusion

FusionInspector.v2.6.0/FusionInspector --fusions $outdir/star-fusion/star-fusion.fusion_predictions.tsv \
                --genome_lib star-fusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir \
                --left_fq $fq1 \
                --right_fq $fq2  \
                -O  $outdir/FusionInspector  \
                --out_prefix finspector \
                --vis
EMBOSS-6.6.0/emboss/transeq  $outdir/FusionInspector/${sample}.finspector.fa  -frame=3 $outdir/${sample}.fusion.fasta
