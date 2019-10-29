import os
import pandas as pd
import sys
import re

phenotype_matrix = sys.argv[1]
outfolder = sys.argv[2]
pheno_to_exclude = sys.argv[3]
covariate_file = sys.argv[4]

phenotag = re.search('([A-Za-z0-9._]+quantif)',phenotype_matrix).group(1)


if not os.path.exists(outfolder):
	os.makedirs(outfolder)
if not os.path.exists(outfolder+'/OUT/'):
	os.makedirs(outfolder+'/OUT/')
if not os.path.exists(outfolder+'/ERR/'):
	os.makedirs(outfolder+'/ERR/')


vcf_file = '/data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells/EGAD00001002663/mergedSGX/BLUEPRINT_SGX_merged.vcf.gz'

#genotype: may need to remove some PCs, to check

#com1 = 'bsub -K -o '+outfolder+'/OUT/genotype_pca.out -e '+outfolder+'/ERR/genotype_pca.err -J genotype_pca -q priority \'~/bin/QTLtools1.1/bin/QTLtools pca --vcf '+vcf_file+' --scale --center --maf 0.05 --distance 50000 --out '+outfolder+'/genotype_pca\''
#print(com1)
#os.system(com1)

#phenotype
com2 = 'bsub -K -o '+outfolder+'/OUT/'+phenotag+'_pca.out -e '+outfolder+'/ERR/'+phenotag+'_pca.err -J '+phenotag+'_pca -q priority \'~/bin/QTLtools1.1/bin/QTLtools pca --bed '+phenotype_matrix+' --scale --center --out '+outfolder+'/'+phenotag+'_pca\''
print(com2)
os.system(com2)

#covariates
covariate_data = pd.read_table(covariate_file,sep='\t',header=None, index_col=0)
covariate_data = covariate_data.transpose()
covariate_data = covariate_data.assign(SampleID="Sex")

#genopca = pd.read_table(outfolder+'/genotype_pca.pca',sep=' ')
phenopca = pd.read_table(outfolder+'/'+phenotag+'_pca.pca',sep=' ')

#sys.exit()

#run qtltools
nperm = 1000
nchunks = 100
covariates = range(-1,50,2)
for i in covariates:
	if i>=0:
		genopheno = phenopca.ix[0:i,:]
		genopheno = pd.concat([genopheno,covariate_data],join="inner")
	else:
		genopheno = covariate_data

	genopheno.to_csv(outfolder+'/covariates_'+str(i+1)+'.txt',sep='\t',float_format='%.3f',index=False)
	for j in range(1,nchunks+1):
		bsubcom1 = 'bsub -o '+outfolder+'/OUT/cov_'+str(i+1)+'_chunk'+str(j)+'.out -e '+outfolder+'/ERR/cov_'+str(i+1)+'_chunk'+str(j)+'.err -J '+phenotag+'_cov_'+str(i+1)+'_chunk'+str(j)+' -q normal -M 16000000 '
		com3 = '\'/home/grey2/bin/QTLtools1.1/bin/QTLtools cis --normal --vcf '+vcf_file+' --bed '+phenotype_matrix+' --cov '+outfolder+'/covariates_'+str(i+1)+'.txt --permute '+str(nperm)+' --chunk '+str(j)+' '+str(nchunks)+' --exclude-phenotypes '+pheno_to_exclude+' --out '+outfolder+'/permutations_cov'+str(i+1)+'_'+str(j)+'_'+str(nchunks)+'.txt\''
		print(bsubcom1+com3)
		os.system(bsubcom1+com3)

#assemble results
for k in covariates:
	com4 = '\'cat '+outfolder+'/permutations_cov'+str(k+1)+'_*_'+str(nchunks)+'.txt | gzip -c > '+outfolder+'/permutations_cov'+str(k+1)+'_full.txt.gz; '
	com5 = '/software/R/3.4.2/bin/Rscript /home/grey2/bin/QTLtools1.1/script/runFDR_cis.R '+outfolder+'/permutations_cov'+str(k+1)+'_full.txt.gz 0.05 '+outfolder+'/permutations_cov'+str(k+1)+'\''
	bsubcom2 = 'bsub -o '+outfolder+'/OUT/cov_'+str(k+1)+'_assemble.out -e '+outfolder+'/ERR/cov_'+str(k+1)+'_assemble.err -J cov_'+str(k+1)+'_assemble -q priority -w \'ended(\''+phenotag+'*\')\' '
	com6 = bsubcom2+com4+com5
	print(com6)
	os.system(com6)


