#!/bin/bash

CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
DATADIR_0=/home/users/a/avalosma/scratch/0_CRD
DATADIR=/home/users/a/avalosma/scratch/1_CRD
OUT_FOLDER=/home/users/a/avalosma/scratch/2_CRD
mkdir -p $OUT_FOLDER $OUT_FOLDER/OUT

#Loop on all cell types:
# declare -a Cell_Types=("EGAD00001002673" "EGAD00001002672" "EGAD00001002670")

cell_type="EGAD00001002670"

for c in $(seq 1 2); do
    echo "$c"
    LI=$DATADIR_0/${cell_type}_merged_residuals.bed.gz
    LM=$DATAFOLDER/${cell_type}_ALL.modules.MOD1.NRE2.txt.gz
    LT=$DATADIR/$cell_type\.chr$c\.module.txt.gz
    LO=$OUT_FOLDER/$cell_type\.chr$c\.pc1

    cmd="$CLOMICS quantify --bed $LI --region $c --tree $LT $LM --out $LO\.pc1.txt.gz --pca 1 --normal"
    cmd2="$CLOMICS quantify --bed $LI --region $c --tree $LT $LM --out $LO\.mean.txt.gz --mean --normal"
    cmd3="$CLOMICS quantify --bed $LI --region $c --tree $LT $LM --out $LO\.loom.txt.gz --loo 0 --normal"
    
    echo "$cmd"
    JOB=$OUT_FOLDER/OUT/${cell_type}_${c}_quantify
    sbatch -J $JOB\.job --partition=mono-EL7 --time=00:08:00 -o $JOB\.out -e $JOB\.err --wrap="$cmd && $cmd2 && $cmd3" 
done
