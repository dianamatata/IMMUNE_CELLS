#!/bin/bash

### Mapping of large sets of high-throughput sequencing reads to a reference genome is one of the foundational steps in RNA-seq data analysis. 
### The STAR software package performs this task with high levels of accuracy and speed.

# we have .bam files that have been processed with STAR.
# file name type: S013DLB3.total_RNA.paired_end_STAR_wtsi.GRCh37.20150724.bam
# it is aligned with ref genome GRCh37, the STAR mapping quality is 255 by default

# G code : /home/grey2/bin/QTLtools1.1/bin/QTLtools quan --gtf ../reference_files/gencode.v15.annotation.gtf.gz --bam ./*.bam --out qtltools_quantification --filter-mapping-quality 255 --rpkm "
# folders with RNA seq data in FOLDERS
# https://qtltools.github.io/qtltools/ to check the quan fct
# https://www.gencodegenes.org/human/release_33lift37.html
# --filter-mismatch-total is 8 for 50bp, 12 for 75bp, 16 for 100bp

BAM_FOLDER=/home/users/a/avalosma/scratch//Blueprint/RNA_seq
OUT_FOLDER=/home/users/a/avalosma/scratch//Blueprint/RNA_seq/quantification
GENCODE_FOLDER=/home/users/a/avalosma/scratch/Blueprint/ref_files

mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT

#for folder in 'EGAD00001002675' 'EGAD00001002671' 'EGAD00001002674' ; do
for folder in '/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells/EGAD00001002675' ; do
	for gtf_file in 'gencode.v15.annotation.gtf.gz' 'gencode.v33lift37.annotation.gtf.gz' ; do
		# cmd="QTLtools quan --bam $BAM_FOLDER/$folder/*.bam --gtf $GENCODE_FOLDER/$gtf_file --out $OUT_FOLDER/qtltools_quantification_${folder}_${gtf_file} --filter-mapping-quality 255 --filter-mismatch-total 16 --rpkm"
                cmd="QTLtools quan --bam $folder/*.bam --gtf $GENCODE_FOLDER/$gtf_file --out $OUT_FOLDER/qtltools_quantification_${folder}_${gtf_file} --filter-mapping-quality 255 --filter-mismatch-total 16 --rpkm"	
		echo "$cmd"
		#eval $cmd
		#FILENAME=qtl_quan_${folder}_${gtf_file}_job
		FILENAME='filename'
		#sbatch -J $FILENAME.job --partition=mono-EL7 --time=02:00:00 -o $OUT_FOLDER/OUT/$FILENAME.out -e $OUT_FOLDER/OUT/$FILENAME.err --wrap="$cmd"
		sbatch -J $FILENAME.job --partition=mono-EL7 --time=02:00:00 --wrap="$cmd"
	done
done
