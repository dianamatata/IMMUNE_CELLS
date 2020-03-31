#!/bin/bash

# goal: compute PCs

BED=/data/unige/funpopgen/davalos/project/GEUVADIS/step1_preparing_VCF_BED/bed_10percent.bed.gz
PCA_FOLDER=/data/unige/funpopgen/davalos/project/GEUVADIS/step2_PCA
PCA_FILE=$PCA_FOLDER/bed_10percent.pca

# apply PCA
echo "performing QTLtools PCA"
QTLtools pca --bed $BED --scale --center --out $PCA_FOLDER/bed_10percent

# create PC covariate files
echo "create PC covariate files"

for PC in 1 2 5 10 20 30 40 50; do
    echo create PCA file $PC
    echo "head -$((PC +1)) $PCA_FILE > $PCA_FOLDER/PCs/PC_$PC.txt"
    head -$((PC +1)) $PCA_FILE > $PCA_FOLDER/PCs/PC_$PC.txt
done

