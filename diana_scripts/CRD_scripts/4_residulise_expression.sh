#!/bin/bash

QTLTOOLS=/home/users/a/avalosma/bin/qtltools
OUTFOLDER=/home/users/a/avalosma/scratch/3_CRD
mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT

# for RNA
cell_type='EGAD00001002675'

for PC in 0 10 20 30 40 50; do
	# to modify
	INP=/data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells/${cell_type}/qtltools_quantification_filtered.gene.rpkm.bed.gz
	COV=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD//THREE_CELL_TYPES/QTL_TOOLS/${cell_type}_v3.0/covariates_$PC.txt
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD//THREE_CELL_TYPES/QTL_TOOLS/${cell_type}_v3.0/${cell_type}_RNA.PC$PC\.bed
	
	cmd="$QTLTOOLS correct --bed $INP --cov $COV --normal --out $OUT && $BGZ $OUT && $TBX -p bed $OUT\.gz"
	JOB=residualize_$PC
	sbatch -J $JOB\.job --partition=mono-EL7 --time=00:01:00 -o $OUT_FOLDER/OUT/$JOB.out -e $OUT_FOLDER/OUT/$JOB.err --wrap="$cmd"
done
