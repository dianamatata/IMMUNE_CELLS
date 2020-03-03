#!/bin/bash

CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
DATADIR_0=/home/users/a/avalosma/scratch/0_CRD
DATADIR=/home/users/a/avalosma/scratch/1_CRD
OUT_FOLDER=/home/users/a/avalosma/scratch/2_CRD
mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
	LI=$DATADIR_0/${cell_type}_merged_residuals.bed.gz
        LM=$DATADIR/${cell_type}_ALL.modules.MOD1.NRE2.txt.gz
	for c in $(seq 1 222); do
		echo "$cell_type $c"
		LT=$DATADIR/$cell_type\.chr$c\.module.txt
		LO=$OUT_FOLDER/$cell_type\.chr$c\

		cmd="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO\.pc1.txt.gz --pca 1 --normal"
		cmd2="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO\.mean.txt.gz --mean --normal"
		cmd3="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO\.loom.txt.gz --loo 0 --normal"

		#echo "$cmd"
		JOB=$OUT_FOLDER/OUT/${cell_type}_chr${c}_quantify
		sbatch -J $JOB\.job --partition=mono-EL7 --time=00:01:00 -o $JOB\.out -e $JOB\.err --wrap="$cmd && $cmd2 && $cmd3" 
	done
done
