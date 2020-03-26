#!/bin/bash

CLOMICs=/home/users/a/avalosma/bin/clomics/bin/clomics
DATADIR=/home/users/a/avalosma/scratch/0_CRD
DATADIR1=/home/users/a/avalosma/scratch/1_CRD
OUPUT_FOLDER=/home/users/a/avalosma/scratch/Gene_CRD
OUTFOLDER=/home/users/a/avalosma/scratch/Gene_CRD/quantify
mkdir -p $OUT_FOLDER $OUTFOLDER $OUTFOLDER/OUT

for cell_type_data in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
	for cell_type_crd in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
		if [[ "$cell_type_data" != "$cell_type_crd" ]] ; then
			name=$(echo ${cell_type_data: (-2)})data_vs_$(echo ${cell_type_crd: (-2)})crd
			echo $name
			LI=$DATADIR/${cell_type_data}_merged_residuals.bed.gz	
			LM=$DATADIR1/${cell_type_crd}_ALL.modules.MOD1.NRE2.txt.gz
	
			for c in $(seq 1 22); do
				LO=$OUTFOLDER/${name}.chr${c}
				LT=$DATADIR1/${cell_type_crd}.chr${c}.module.txt.gz
			
		                cmd="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.pc1.txt.gz --pca 1 --normal"
                		cmd2="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.mean.txt.gz --mean --normal"
               			cmd3="$CLOMICs quantify --bed $LI --region $c --tree $LT $LM --out $LO.loom.txt.gz --loo 0 --normal"
                		echo $cmd; eval $cmd; eval $cmd2; eval $cmd3
				#sbatch -J ${name}_chr${c}_quantify\.job --partition=mono-EL7 --time=00:01:00 --wrap="$cmd && $cmd2 && $cmd3"

			done
                fi
        done
done
