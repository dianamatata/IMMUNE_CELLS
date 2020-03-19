#!/bin/bash
# goal: create .peak file

CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
DATADIR=/home/users/a/avalosma/scratch/1_CRD
OUT_FOLDER=/home/users/a/avalosma/scratch/3_CRD
mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
	LM=$DATADIR/${cell_type}_94s_ALL.modules.MOD1.NRE2.txt.gz
	for c in $(seq 1 2); do
                LT=$DATADIR/${cell_type}_94s.chr$c\.tree.txt.gz
                LO=$OUT_FOLDER/${cell_type}_94s.chr$c\.peak
		cmd="$CLOMICs gpeak --tree $LT --mod $LM --out $LO\.txt.gz"
		JOB=${cell_type}_chr${c}_peak
		sbatch -J $JOB\.job --partition=mono-EL7 --time=00:01:00 --wrap="$cmd"
	done
done

