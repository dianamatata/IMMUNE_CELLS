#!/bin/bash

# Goal: extract the 7th column of each file and stack them together
# to obtain https://qtltools.github.io/qtltools/pages/mode_quant.html , the last table they display in this page

DATA_FOLDER=/home/users/a/avalosma/scratch/RNA_quantify_indep

for cell_type in 'EGAD00001002671'  'EGAD00001002674' 'EGAD00001002675' ; do
	f=$(echo "file_${cell_type}.txt")
	
	# get the 7th line of all the files belonging to the same cell type and paste in line (-s)
	for i in ${cell_type}*.gene.rpkm.bed; do cat $i | cut -f7 | paste -s; done > $f
	
	# transpose file
	cat $f | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' > temp.txt
	
	# take the sample name of the files, and change the header line
	cat temp.txt | head -1 | tr '\t' '\n' | rev | cut -d '/' -f1 | rev | cut -d '_' -f1 | tr '\n' '\t' > temp2.txt
	echo "" >> temp2.txt 
	sed '1d' temp.txt >> temp2.txt
	cut -c 2- temp2.txt > $f
	rm temp.txt temp2.txt
done


#2nd method
for cell_type in 'EGAD00001002671'  'EGAD00001002674' 'EGAD00001002675' ; do
	paste ${cell_type}*.gene.rpkm.bed | awk '{printf "%s",$1; for (i=7;i<=NF;i=i+7) printf " %s",$i; print ""}' > file1
	#to do first line short names
	#cat file1 | head -1 | tr ' ' '\n' | rev | cut -d '/' -f1 | rev | cut -d '_' -f1 |  tr '\n' '\t'
	paste <(cat EGAD00001002674_qltools_quan_S00FK4B5.GnJ6geaB7TH.gene.rpkm.bed | cut -f1-6) file1 > file2
done


