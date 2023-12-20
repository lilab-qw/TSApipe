#!/bin/bash

# User-defined parameters

PROJECT_FILE=$1 #/u/laumontc/tsaPaper/${SAMPLE_NAME}/kmerFilt/${SAMPLE_NAME}.${Lcount}.mtechumanAll.0.project
KMER_FILE=$2 #/u/laumontc/tsaPaper/${SAMPLE_NAME}/kmerFilt/${SAMPLE_NAME}.${Lcount}.mtechumanAll.0.Lcount
K_LEN=33
ASS_TYPE='linear' # graph --> defines the type of assembly you wanna do
CHECK_RC='FALSE' # TRUE --> check the reverse complement during assembly or not
PATH_OUT=$3
Lcount=$4

# Create relevant arborescence to store results
ASS_NAME=`echo "${KMER_FILE%.Lcount}" | awk -F [/] '{print $NF}'`

if [ ! -d "$PATH_OUT/" ]; then
    # echo 'Does not exists'
    mkdir $PATH_OUT
    cd $PATH_OUT
    mkdir $ASS_TYPE'_'$ASS_NAME
    cd $ASS_TYPE'_'$ASS_NAME
else
    # echo 'Exists'
    cd $PATH_OUT
    mkdir $ASS_TYPE'_'$ASS_NAME
    cd $ASS_TYPE'_'$ASS_NAME
fi


# only keep kmer and cancer Lcount
# assembly does not like to have too many columns
cut -f2-3 $KMER_FILE > ./$ASS_NAME.assembly


if [ $ASS_TYPE == "linear" ] && [ $CHECK_RC == 'FALSE' ]; then
    echo 'Starting linear assembly - No reverse complement checking';
    kmer_assembly -p $PROJECT_FILE -k ./$ASS_NAME.assembly \
                                                -o $PATH_OUT/$ASS_TYPE'_'$ASS_NAME/  -l $K_LEN;

elif [ $ASS_TYPE == "linear" ] && [ $CHECK_RC == 'TRUE' ]; then
    echo 'Starting linear assembly - Reverse complement checking';
    kmer_assembly -r -p $$PROJECT_FILE -k ./$ASS_NAME.assembly \
                                                -o $PATH_OUT/$ASS_TYPE'_'$ASS_NAME/  -l $K_LEN;

elif [ $ASS_TYPE == "shortestPath" ] && [ $CHECK_RC == 'FALSE' ]; then
    echo 'Starting graph assembly - No reverse complement checking';
    kmer_assembly -p $$PROJECT_FILE -k ./$ASS_NAME.assembly \
                                                -o $PATH_OUT/$ASS_TYPE'_'$ASS_NAME/  -l $K_LEN;

elif [ $ASS_TYPE == "shortestPath" ] && [ $CHECK_RC == 'TRUE' ]; then
    echo 'Starting graph assembly - Reverse complement checking';
    kmer_assembly -r -p $$PROJECT_FILE -k ./$ASS_NAME.assembly \
                                                -o $PATH_OUT/$ASS_TYPE'_'$ASS_NAME/  -l $K_LEN;

else
    echo 'ASS_TYPE variable must be set to "linear" or "shortestPath"';
    echo 'CHECK_RC variable must be set to "TRUE" or "FALSE"';
fi
