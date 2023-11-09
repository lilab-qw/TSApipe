sample=$1
fq1=$2
fq2=$3
outdir=$4
NB_THREADS=4
jellyfish=/home/wangmengyao/miniconda2/bin/jellyfish
fastx_reverse_complement=/home/qinyurong/biosoft/bin/fastx_reverse_complement

# jellyfish parameters
STRANDED='FALSE' # or TRUE
LEN=(24 33)
HASH_SIZE='1G'
set +o posix
if [[ $STRANDED = 'TRUE' ]]; then
     # Reverse Complement R1 and stranded jf database
     echo `date "+%Y-%m-%d %H:%M:%S":` 'Branch for stranded data...'
     echo '    Reverse complementing R1 fastq...'  `date "+%Y-%m-%d %H:%M:%S"`
     zcat $fq1 | $fastx_reverse_complement -z $outdir/$sample\_R1.rc.fastq.gz

     echo '    Generating non-canonical jf database and other info...' `date "+%Y-%m-%d %H:%M:%S"`
     for l in "${LEN[@]}"
     do
        echo '        Length:' $l'...' `date "+%Y-%m-%d %H:%M:%S"`

        $jellyfish count -m $l -s $HASH_SIZE -t $NB_THREADS -o $outdir/$sample.trim.R1rc.$l.jf <  zcat $outdir/$sample\_R1.rc.fastq.gz $fq2

        $jellyfish histo $outdir/$sample.trim.R1rc.$l.jf > $outdir/$sample.trim.R1rc.$l.histo
        $jellyfish info $outdir/$sample.trim.R1rc.$l.jf > $outdir/$sample.trim.R1rc.$l.info
     done
else

     # No reverse complement and unstranded jf database
     echo `date "+%Y-%m-%d %H:%M:%S":` 'Branch for unstranded data...'
     echo '    Generating canonical jf database (-C) and other info...' `date "+%Y-%m-%d %H:%M:%S"`
     for l in "${LEN[@]}"
     do
         echo '        Length:' $l'...' `date "+%Y-%m-%d %H:%M:%S"`
         $jellyfish count -m $l -s $HASH_SIZE -t $NB_THREADS -C -o $outdir/$sample.trim.$l.jf <(zcat $fq1) <(zcat $fq2)
         $jellyfish histo $outdir/$sample.trim.$l.jf > $outdir/$sample.trim.$l.histo
         $jellyfish info $outdir/$sample.trim.$l.jf > $outdir/$sample.trim.$l.info
 done
fi
