#!/bin/bash

#pull read names from every .fa file in the directory 
DIR=/path/to/directory/
#lists the path to all files we are going to use
for f in ${DIR}*.fq
do
grep -o '^@[^ ]*\b' $f >> read_IDs_1.txt
done

cut -c 2- read_IDs_1.txt > read_IDs.txt
