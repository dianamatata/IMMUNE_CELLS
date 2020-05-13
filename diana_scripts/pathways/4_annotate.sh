#!/bin/bash

$inFolder='trans_hubs_genes'
$refGen='hg37'
for file in "$inFolder"/*; do
        printf $file
	file_nbr=$(echo $file | awk -F '[_.]'  '{print $4"_"$5}')
	printf $file_nbr
	annotatePeaks.pl $file $refGen > ${outFile}_$file_nbr
done

