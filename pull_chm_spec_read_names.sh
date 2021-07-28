#!/bin/bash

#pull read names from every .fa file in the directory 
#make one file of read names for each .fa file
DIR=/path/to/directory/
#lists the path to all files we are going to use
for f in ${DIR}*.fq
do
var=$(echo "${f##*/}" | cut -f1 -d".")
grep -o '^@[^ ]*\b' $f >> chm_spec_reads/${var}_IDs_1.txt
cut -c 2- chm_spec_reads/${var}_IDs_1.txt > chm_spec_reads/${var}_IDs.txt
echo $var
done


