#!/bin/bash

DATADIR=/home/users/a/avalosma/scratch/2_CRD

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do

	# unzip all the mean.txt.gz files
	gunzip $DATADIR/${cell_type}.chr*.mean.txt.gz

	# take the first file fully and then everything except the first line into unsorted 
	awk 'FNR>1||NR==1' $DATADIR/${cell_type}.chr*.mean.txt > $DATADIR/${cell_type}.ALLchr.mean.unsorted.txt
	#sort all the lines except header, by col 1 first and then col 2 by n numerical value
	awk 'NR == 1; NR > 1 {print $0 | "sort -V -k1,1 -k2,2n"}' $DATADIR/${cell_type}.ALLchr.mean.unsorted.txt | bgzip -c > $DATADIR/${cell_type}.ALLchr.mean.txt.gz

	tabix -p bed $DATADIR/${cell_type}.ALLchr.mean.txt.gz
	bgzip $DATADIR/${cell_type}.chr*.mean.txt
done
