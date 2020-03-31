#!/bin/bash

# Goal: compute mQTLs

DATAFOLDER=/home/users/a/avalosma/scratch/Blueprint/DNA_meth/EGAD00010000850/processed_data
OUTFOLDER=/home/users/a/avalosma/scratch/mQTLs
mkdir -p $OUTFOLDER $OUTFOLDER/OUT $OUTFOLDER/ERR

phenotype_matrix=$DATAFOLDER/mono_meth_M_20151028_formatted.bed.gz
#pheno_to_exclude=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/methylomics/mono_phenotypes_to_exclude.txt #is empty
covariate_file=/home/users/a/avalosma/scratch/Blueprint/ref_files/Donor_id_Sex_merged.map
vcf_file=/home/users/a/avalosma/scratch/Blueprint/ref_files/EGAD00001002663/All_chr.BPWP10_13_12_15.vcf.gz
# where does this vcf comes from? is it in the blueprint ?

cmd_old='bsub -K -o '+outfolder+'/OUT/genotype_pca.out -e '+outfolder+'/ERR/genotype_pca.err -J genotype_pca -q priority \'~/bin/QTLtools1.1/bin/QTLtools pca --vcf '+vcf_file+' --scale --center --maf 0.05 --distance 50000 --out '+outfolder+'/genotype_pca\'

cmd_old_2=com2 = 'bsub -K -o '+outfolder+'/OUT/'+phenotag+'_pca.out -e '+outfolder+'/ERR/'+phenotag+'_pca.err -J '+phenotag+'_pca -q priority \'~/bin/QTLtools1.1/bin/QTLtools pca --bed '+phenotype_matrix+' --scale --center --out '+outfolder+'/'+phenotag+'_pca\'


cmd=""
#eval $cmd

cmd2=""
#eval $cmd2

nperm=1000
nchunks=100

for i in $(seq -1 2 50); do
	if [ $i -ge 0 ];
	then
	        print $i
	else
	fi 
done

