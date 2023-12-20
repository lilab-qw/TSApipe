sample=$1 
mutation_type=$2
MS_filter_file=$3
hla=$4
outdir=$5
cd $outdir
awk -F '\t' '$4!="Sq"   {print$4}' ${outdir}/$MS_filter_file  > ${outdir}/${sample}.${mutation_type}.MAPS-1.fasta
sort  ${sample}.${mutation_type}.MAPS-1.fasta |uniq  > ${sample}.${mutation_type}.MAPS-2.fasta
awk -v x=7 '{if(length($1)>x){print $0}}' ${sample}.${mutation_type}.MAPS-2.fasta > ${sample}.${mutation_type}.MAPS.fasta
rm ${sample}.${mutation_type}.MAPS-1.fasta ${sample}.${mutation_type}.MAPS-2.fasta 
python 04.Identification/length.py ${outdir}/$MS_filter_file  ${outdir}/${sample}.${mutation_type}-1.csv
cat   ${outdir}/${sample}.${mutation_type}-1.csv |awk '!a[$4]++{print}'> ${sample}.${mutation_type}.csv
rm ${sample}.${mutation_type}-1.csv
echo $(date) ${sample}.${mutation_type}.csv finished ...

for hla in `cat $hlaI`
do
        netMHCpan-4.0/netMHCpan -a $hla -p -inptype 1 -l 8,9,10,11 -s -f ${sample}.${mutation_type}.MAPS.fasta  -BA > $outdir/${sample}.$hla.binding.xls
done
echo $(date) ${sample}.MHCI.$hla.binding.xls  finished ...

ls $outdir/${sample}.*.binding.xls > $outdir/${sample}.MHCI.binding.list

perl 04.Identification/combine.MHC.binding.pl $outdir/${sample}.MHCI.binding.list $outdir/${sample}.MHCI.binding.merge.txt
echo $(date) ${sample}.MHCI.binding.merge.txt  finished ...

perl 04.Identification/combineMS.MHC.pl ${outdir}/${sample}.${mutation_type}.csv   $outdir/${sample}.MHCI.binding.merge.txt  $outdir/$sample.capaMHC.MHCbinding.csv
awk -F '\t' '{if( $(NF-1) == "Min_rank" ||$(NF-1)<2) {print}}' ${sample}.capaMHC.MHCbinding.csv > ${sample}.capaMHC.MHCbinding.final.csv
sed -i 's/\r//g' ${sample}.capaMHC.MHCbinding.final.csv
echo $(date) $sample.capaMHC.MHCbinding.csv finished ...
