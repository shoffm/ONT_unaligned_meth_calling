#!/bin/bash

#load fast5 api
module load ont-fast5-api

#input path
OUTPATH_SUB=/data/hoffmansf/read_isolation/fast5_subset_output
READ_ID_LIST=/home/hoffmansf/read_isolation/read_IDs.txt

#cat /data/hoffmansf/read_isolation/fast5_paths_2_lin.txt | while read line
cat /data/hoffmansf/read_isolation/fast5_paths_4_lin.txt | while read line
do
#echo $line
#NAME=$(echo $line | awk -F'/' ' { print $NF } ')
NAME=$(echo "${line//\//.}")
NAME2="${NAME:1}"
#echo ${NAME2}
#ID=${$line | awk -F'/' '{ print $NF } '}
OUTPATH=${OUTPATH_SUB}/${NAME2}/
#echo ${OUTPATH}
#INPATH=${line%/*}
#echo $line
fast5_subset --input ${line} --save_path ${OUTPATH} --read_id_list ${READ_ID_LIST}
done

