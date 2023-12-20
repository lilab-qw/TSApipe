#!/bin/bash

# Sample info used for query
SAMPLE_NAME=$1
PATH_TO_JF=$2

# Minimal count for cancer k-mers (>=)
LCOUNT=$3 #7
# Normal k-mer database(s) that should be used for filtering
PATH_TO_JF_FILTER=$4
# Path to write output file
PATH_OUT=$5

echo `date "+%Y-%m-%d %H:%M:%S":` 'Dumping cancer k-mers having an occurrence superior or equal to '$LCOUNT' in sample '$PATH_TO_JF'...'
jellyfish dump -c -t -L $LCOUNT -o $PATH_OUT/$SAMPLE_NAME.$LCOUNT.tmp $PATH_TO_JF
# add ids for each k-mers - useful in the worflow
awk '{print NR "\t" $0}' $PATH_OUT/$SAMPLE_NAME.$LCOUNT.tmp > $PATH_OUT/$SAMPLE_NAME.$LCOUNT.txt
# remove temporary file (that does not contain id for each k-mer)
rm -v $PATH_OUT/$SAMPLE_NAME.$LCOUNT.tmp
echo `date "+%Y-%m-%d %H:%M:%S":` 'Done!'
echo ''


echo `date "+%Y-%m-%d %H:%M:%S":` 'Querying k-mers in normal k-mer database(s)...'
# cut -f2 $PATH_OUT/$SAMPLE_NAME.$LCOUNT.txt > $PATH_OUT/queries
count=1
for jf in "${PATH_TO_JF_FILTER[@]}"
do
    echo '    querying in: '$jf'...'
    cut -f2 $PATH_OUT/$SAMPLE_NAME.$LCOUNT.txt | jellyfish query -i $jf > $PATH_OUT/$count.tmp
    # if [ "$count" -eq "${#PATH_TO_JF_FILTER[@]}" ]; then
    #   header+=${FILTER_NAME[$count-1]}
    # else
    #   header+=${FILTER_NAME[$count-1]}'\t'
    # fi
    count=$((count+1))
done
echo `date "+%Y-%m-%d %H:%M:%S":` 'Done!'
echo ''



echo `date "+%Y-%m-%d %H:%M:%S":` 'Generating the allcount file...'
# works for as many file as you want ! sort -V fixes the 10.tmp 1.tmp order pb!
paste $PATH_OUT/$SAMPLE_NAME.$LCOUNT.txt $(find $PATH_OUT -name "*.tmp" | sort -V) > $PATH_OUT/$SAMPLE_NAME.$LCOUNT.allcount
# remove temporary files (since counts are all concatenated in the .allcount)
rm -v $PATH_OUT/*.tmp
# rm queries
rm -v $PATH_OUT/$SAMPLE_NAME.$LCOUNT.txt
echo `date "+%Y-%m-%d %H:%M:%S":` 'Done!'
echo ''



echo 'Filtering k-mers...'
# echo $(wc -l $PATH_OUT/test2) > $PATH_OUT/$SAMPLE_NAME$TYPE$LEN.$LCOUNT.filtering
# count the number of initial cancer k-mers and save it in .filtering file
echo $(wc -l $PATH_OUT/$SAMPLE_NAME.$LCOUNT.allcount) > $PATH_OUT/$SAMPLE_NAME.$LCOUNT.filtering

echo $(wc -l $fileName.$colName.0) >> $PATH_OUT/$SAMPLE_NAME$TYPE$LEN.$LCOUNT.filtering


# removing FILTER_NAME
colToFilt=''
colName=''
fileName=''
header="id\tkmer\t${SAMPLE_NAME}"

for i in "${!PATH_TO_JF_FILTER[@]}"
do
    colToFilt=$((i + 1 + 3)) # +1 to correct for 0-based, +3 because files already contains 3 fixed columns (id, k-mers, cancerCount)
    tmpColName=`echo "${PATH_TO_JF_FILTER[$i]}" | awk -F [/] '{print $NF}'`
    colName=`echo "$tmpColName" | awk -F [.] '{print $1}'`

    header+='\t'${colName}
    # echo $header
    # echo $colToFilt $colName
    if [ "$colToFilt" -eq 4 ]; then
        # awk -v c=$colToFilt '$c == 0 {print $0}' $PATH_OUT/test2 > $PATH_OUT/$SAMPLE_NAME$TYPE$LEN.$LCOUNT.$colName.0
        awk -v c=$colToFilt '$c == 0 {print $0}' $PATH_OUT/$SAMPLE_NAME.$LCOUNT.allcount > $PATH_OUT/$SAMPLE_NAME.$LCOUNT.$colName.0
        echo $(wc -l $PATH_OUT/$SAMPLE_NAME.$LCOUNT.$colName.0) >> $PATH_OUT/$SAMPLE_NAME.$LCOUNT.filtering
    else
        tmpColNameP=`echo "${PATH_TO_JF_FILTER[$i-1]}" | awk -F [/] '{print $NF}'`
        colNameP=`echo "$tmpColNameP" | awk -F [.] '{print $1}'`
        fileName=$(find $PATH_OUT -name "$SAMPLE_NAME\.$LCOUNT.*$colNameP\.0")
        awk -v c=$colToFilt '$c == 0 {print $0}' $fileName > $fileName.$colName.0
        echo $(wc -l $fileName.$colName.0) >> $PATH_OUT/$SAMPLE_NAME.$LCOUNT.filtering
    fi
done


# Add header to all files
echo `date "+%Y-%m-%d %H:%M:%S":` 'Adding header to file...!'
sed -i "1i ${header}" $PATH_OUT/$SAMPLE_NAME.$LCOUNT.allcount
# sed -i "1i ${header}" $PATH_OUT/test2

for f in $(find $PATH_OUT -name "$SAMPLE_NAME\.$LCOUNT.*\.0")
do
    sed -i "1i ${header}" $f
    mv $f $f.count
done

# Final cleanup
echo 'Removing filter jellyfish to RamDisk'
rm $PATH_TO_JF_FILTER
echo `date "+%Y-%m-%d %H:%M:%S":` 'Done!'
echo ''
