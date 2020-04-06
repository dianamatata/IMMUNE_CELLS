#!/bin/bash

# Goal: compute mQTLs

DATAFOLDER=/home/users/a/avalosma/scratch/Blueprint/DNA_meth/EGAD00010000850/processed_data
OUTFOLDER=/home/users/a/avalosma/scratch/mQTLs
mkdir -p $OUTFOLDER $OUTFOLDER/OUT $OUTFOLDER/ERR

## comments
#pheno_to_exclude=/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/methylomics/mono_phenotypes_to_exclude.txt
#is empty
# genotype: NO NEED to remove PCs, no structure observed
# where does this vcf comes from? is it in the blueprint ?


# compute PCA wth QTLtools * no samples to exclude * MAF greater than 0.05 * Sites every 50000 bp #samples = 197 
# Running time: 569 sec

covariate_file=/home/users/a/avalosma/scratch/Blueprint/ref_files/Donor_id_Sex_merged.map
vcf_file=/home/users/a/avalosma/scratch/Blueprint/ref_files/EGAD00001002663/All_chr.BPWP10_13_12_15.vcf.gz
cmd="QTLtools pca --vcf $vcf_file --scale --center --maf 0.05 --distance 50000 --out $OUTFOLDER/genotype_pca"
eval $cmd


for cell_type in 'monocytes' 'neutrophils' 'tcells' ; do 
	phenotype_matrix=$DATAFOLDER/methyldata_${cell_type}.bed.gz
	cmd2="QTLtools pca --bed $phenotype_matrix --scale --center --out $OUTFOLDER/${cell_type}_pca"
	eval $cmd2
done

# transpose part with covariates
cat $covariate_file | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}' | tr ':' '\t' | cut -c 2- > $OUTFOLDER/covariate_file_sexID.txt

phenopca=$OUTFOLDER/${cell_type}_pca


# run qtltools

nperm=1000
nchunks=100

for i in $(seq -1 2 50); do # loop in covariates
	if [ $i -ge 0 ];
	then
	        print $i
	else
	fi 
done

#          sbatch -J ${cell_type}_pca_methyl.jopb --partition=mono-shared-EL7 --time=00:30:00 --wrap="$cmd2"
mkdir -o $OUTFOLDER/covariates $OUTFOLDER/covariates/OUT $OUTFOLDER/covariates/ERR

for j in $(seq 1 1 $nchunks+1); do # chunks
	cmd3="QTLtools cis --normal --vcf $vcf_file --bed $phenotype_matrix --cov $OUTFOLDER/covariates_$(i+1).txt --permute $nperm  --chunk $j $nchunks --out $OUTFOLDER/covariates/permutations_cov_$(i+1)_${j}_${nchunks}.txt"

	eval $cmd3
done

# assemble results
for i in $(seq -1 2 50); do # loop in covariates

com4 = "cat  $OUTFOLDER/covariates/permutations_cov_$(i+1)_${j}_${nchunks}.txt | gzip -c > _full.txt.gz"
com5 = 'Rscript /home/grey2/bin/QTLtools1.1/script/runFDR_cis.R ' + outfolder + '/permutations_cov' + str(
        k + 1) + '_full.txt.gz 0.05 ' + outfolder + '/permutations_cov' + str(k + 1) + '\''
    bsubcom2 = 'bsub -o ' + outfolder + '/OUT/cov_' + str(k + 1) + '_assemble.out -e ' + outfolder + '/ERR/cov_' + str(
        k + 1) + '_assemble.err -J cov_' + str(k + 1) + '_assemble -q priority -w \'ended(\'' + phenotag + '*\')\' '
    com6 = bsubcom2 + com4 + com5'

done

