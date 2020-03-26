#!/bin/bash

# RNA to be updated
RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002675_v3.0/EGAD00001002675_RNA.PC10.bed.gz
MOD=/home/users/a/avalosma/scratch/Gene_CRD/quantify
OUT=/home/users/a/avalosma/scratch/Gene_CRD/mapping_aCRD_gene
mkdir -p $OUT

K=100
for cell in '70' '72' '73' ; do
	for cell2 in '70' '72' '73' ; do
		for k in $(seq 1 $K); do
			name=${cell}data_vs_${cell2}crd
			MOD_NAME=$MOD/${name}.ALLchr.mean.txt.gz
			cmd="QTLtools cis --vcf $MOD_NAME --bed $RNA --permute 200 --chunk $k $K --out $OUT/${name}_gene_CRD_mean_chunk$k\.txt"
			sbatch -J qtltoolscis.job --partition=mono-EL7 --time=01:00:00 --wrap="$cmd"
		done
	done
done

