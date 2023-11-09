#!/bin/bash

set -e

stardir=/home/BIOINFO_TOOLS/alignment_tools/STAR/STAR_020201/bin/Linux_x86_64/
samtools=/home/BIOINFO_TOOLS/mutation_tools/SamTools/SamTools-1.9/samtools
genome=/home/BIOINFO_DATABASE/reference/genome_DNA/Homo_sapiens/hg38/STAR2_AddAnno_index/
annotation=/mnt/delta_DT/TSA/cell_line/RNA-seq/RNA-seq/script/00.Alignment/homo_sapiens/Homo_sapiens.GRCh38.96.chr.gtf

if [ -z "$4" ]
then
    echo "Usage: $0 <sample> <read1> <read2> <workdir>"
    exit 1
fi
sample=$1
r1=$2
r2=$3
dir=$4

### establish local work dir
workdir=$dir/align
mkdir -p ${workdir}

### ALIGNMENT ###
#################
echo "$(date) -- Starting alignment of $sample"

outfile_jct="${dir}/${sample}.junctions"
outfile_aln="${dir}/${sample}.bam"
cd ${workdir}
$stardir/STAR --genomeDir $genome --readFilesIn $r1 $r2 --runThreadN 4 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory  --readFilesCommand zcat --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --sjdbGTFfile $annotation --limitSjdbInsertNsj 2000000 --outSAMunmapped None --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --outSAMattrRGline ID:${sample} SM:${sample} --twopassMode Basic --outSAMmultNmax 1 --limitBAMsortRAM 85420942210  > ${sample}.align.log && \
#--limitBAMsortRAM 70000000000

mv ${workdir}/Aligned.sortedByCoord.out.bam $outfile_aln && \
#mv ${workdir}/SJ.out.tab $outfile_jct && \
$samtools index $outfile_aln && \
echo "$(date) -- Alignment done"

echo "$(date) -- DONE"
