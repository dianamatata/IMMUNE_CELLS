#!/bin/bash

inFolder='/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_genes'
refGen='hg19'
outFile='/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/pathways/annotations/annotate'
for file in "$inFolder"/*; do
        printf $file
	file_nbr=$(echo $file | awk -F '[_.]'  '{print $4"_"$5}')
	printf $file_nbr
	annotatePeaks.pl $file $refGen > ${outFile}_${refGen}_${file_nbr}.txt
done

