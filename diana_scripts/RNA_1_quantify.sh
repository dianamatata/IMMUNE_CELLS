#!/bin/bash

# Goal: This scripts takes the RNA bam files, and call QTLtools quan

DATA_FOLDER=/home/users/a/avalosma/scratch/Blueprint/RNA_seq
REF_FOLDER=/home/users/a/avalosma/scratch/Blueprint/ref_files
OUT_FOLDER=/home/users/a/avalosma/scratch/Blueprint/RNA_quantify
mkdir -p $OUT_FOLDER

# 75 neut, 71 tcell, 74 mono
for cell_type in 'EGAD00001002671'  'EGAD00001002674' 'EGAD00001002675' ; do
	cmd="QTLtools quan --gtf $REF_FOLDER/gencode.v15.annotation.gtf.gz --bam $DATA_FOLDER/$cell_type/*.bam --out $OUT_FOLDER/qltools_quantification --filter-mapping-quality 255 --rpkm"
	echo $cmd
	sbatch -J ${cell_type}_RNAquan.job --partition=mono-EL7 --time=10:00:00 --wrap="$cmd"
	# eval $cmd
done
