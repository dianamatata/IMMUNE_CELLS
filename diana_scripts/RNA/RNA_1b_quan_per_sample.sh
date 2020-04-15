#!/bin/bash

# Goal: This scripts takes the RNA bam files, and call QTLtools quan
# To quantify all samples for which you have RNA-seq, proceed as follows:
#A. Run QTLtools quan for each BAM file separately (i.e. in indepenedent jobs)
#B. Extract the annotation columns using cut -f1-6 on one of the file
#C. Extract the per sample quantifications using cut -f7 for each file separately
#D. Build the entire quantification matrix using the paste to merge B with all quantifications of C
#E. Filter all poorly quantified exons/genes using awk for example


DATA_FOLDER=/home/users/a/avalosma/scratch/Blueprint/RNA_seq
REF_FOLDER=/home/users/a/avalosma/scratch/Blueprint/ref_files
OUT_FOLDER=/home/users/a/avalosma/scratch/RNA_quantify_indep
mkdir -p $OUT_FOLDER

# 75 neut, 71 tcell, 74 mono
for cell_type in 'EGAD00001002671'  'EGAD00001002674' 'EGAD00001002675' ; do
#for cell_type in 'EGAD00001002674'; do
	for file in $DATA_FOLDER/${cell_type}/*.bam ; do
		out_name="$(echo $file | cut -d '/' -f10 | cut -c1-8)"
		cmd="QTLtools quan --gtf $REF_FOLDER/gencode.v15.annotation.gtf.gz --bam $file --out $OUT_FOLDER/${cell_type}_qltools_quan_$out_name --filter-mapping-quality 255 --rpkm"
	        #sbatch -J ${cell_type}_RNAquan.job --partition=mono-EL7 --time=03:00:00 -o $OUTPUT_FOLDER/OUT/${cell_type}_RNAquan_$out_name.out -e $OUTPUT_FOLDER/OUT/${cell_type}_RNAquan_$out_name.err --wrap="$cmd"
		sbatch -J ${cell_type}_RNAquan.job --partition=mono-EL7 --time=03:00:00 --wrap="$cmd"
	done
done
