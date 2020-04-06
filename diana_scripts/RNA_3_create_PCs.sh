#!/bin/bash

# goal: compute PCs

BED_FOLDER=/home/users/a/avalosma/scratch/RNA_rpkm
OUT_FOLDER=/home/users/a/avalosma/scratch/RNA_PC
mkdir -p $OUT_FOLDER

# 75 neut, 71 tcell, 74 mono
for cell_type in 'EGAD00001002671'  'EGAD00001002674' 'EGAD00001002675' ; do
	
	BED=${cell_type}_quantification_filtered.gene.rpkm.bed
        BEDsorted=${cell_type}_quantification_filtered.sorted.gene.rpkm.bed
	# create index file first
	cat $BED_FOLDER/$BED | bedtools sort -header -i > $BED_FOLDER/$BEDsorted
	bgzip $BED_FOLDER/$BEDsorted
	tabix -p bed $BED_FOLDER/${BEDsorted}.gz
	
	# apply PCA
	echo "performing QTLtools PCA"
	QTLtools pca --bed $BED_FOLDER/${BEDsorted}.gz --scale --center --out $OUT_FOLDER/${cell_type}_PCA

	# create PC covariate files
	echo "create PC covariate files"

	for PC in 1 2 5 10 20 30 40 50; do
	     head -$((PC +1)) $OUT_FOLDER/${cell_type}_PCA.pca > $OUT_FOLDER/${cell_type}_PC_$PC.txt
	done
done
