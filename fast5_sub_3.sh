#!/bin/bash

#load fast5 api
module load ont-fast5-api

#input path
OUTPATH_SUB=/path/to/fast5_subset_output
READ_ID_LIST=/path/to/read_IDs.txt

cat /path/to/fast5_paths.txt | while read line
do
NAME=$(echo "${line//\//.}")
NAME2="${NAME:1}"
OUTPATH=${OUTPATH_SUB}/${NAME2}/
fast5_subset --input ${line} --save_path ${OUTPATH} --read_id_list ${READ_ID_LIST}
done

