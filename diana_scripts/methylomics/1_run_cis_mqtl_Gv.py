import os
import pandas as pd
import sys
import re

phenotype_matrix = sys.argv[1]
outfolder = sys.argv[2]
pheno_to_exclude = sys.argv[3]
covariate_file = sys.argv[4]

phenotype_matrix = '/data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells/EGAD00010000850/processed_data/mono_meth_M_20151028_formatted.bed.gz'
outfolder = 'methylomics/QTL_TOOLS/monocytes'
pheno_to_exclude = 'methylomics/mono_phenotypes_to_exclude.txt'
covariate_file = '/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/methylomics/files_test/Donor_id_Sex_merged.map'

phenotag = re.search('([A-Za-z0-9._]+meth)', phenotype_matrix).group(1)  # = 'mono_meth'

vcf_file = '/data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells/' \
           'EGAD00001002663/All_chr.BPWP10_13_12_15.vcf.gz'


# covariates
covariate_data = pd.read_table(covariate_file, sep='\t', header=None, index_col=0)
covariate_data = covariate_data.transpose()
covariate_data = covariate_data.assign(SampleID="Sex")

phenopca = pd.read_table('/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/methylomics/files_test/monocytes_pca.pca', sep=' ')

# sys.exit()

# run qtltools
nperm = 1000
nchunks = 100
covariates = range(-1, 50, 2)
for i in covariates:
    print(i)
    if i >= 0:
        genopheno = phenopca.ix[0:i, :]
        genopheno = pd.concat([genopheno, covariate_data], join="inner")
    else:
        genopheno = covariate_data

    genopheno.to_csv(outfolder + '/covariates_' + str(i + 1) + '.txt', sep='\t', float_format='%.3f', index=False)
    for j in range(1, nchunks + 1):
        bsubcom1 = 'bsub -o ' + outfolder + '/OUT/cov_' + str(i + 1) + '_chunk' + str(
            j) + '.out -e ' + outfolder + '/ERR/cov_' + str(i + 1) + '_chunk' + str(
            j) + '.err -J ' + phenotag + '_cov_' + str(i + 1) + '_chunk' + str(j) + ' -q normal -M 16000000 '
        com3 = '\'/home/grey2/bin/QTLtools1.1/bin/QTLtools cis --normal --vcf ' + vcf_file + ' --bed ' + \
               phenotype_matrix + ' --cov ' + outfolder + '/covariates_' + str(
            i + 1) + '.txt --permute ' + str(nperm) + ' --chunk ' + str(j) + ' ' + str(
            nchunks) + ' --exclude-phenotypes ' + pheno_to_exclude + ' --out ' + outfolder + '/permutations_cov' + str(
            i + 1) + '_' + str(j) + '_' + str(nchunks) + '.txt\''
        print(bsubcom1 + com3)
        os.system(bsubcom1 + com3)

# out:
# bsub -o methylomics/QTL_TOOLS/monocytes/OUT/cov_2_chunk100.out
# -e methylomics/QTL_TOOLS/monocytes/ERR/cov_2_chunk100.err -J mono_meth_cov_2_chunk100 -q normal -M 16000000
# '/home/grey2/bin/QTLtools1.1/bin/QTLtools cis --normal
# --vcf /data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells/EGAD00001002663/All_chr.BPWP10_13_12_15.vcf.gz
# --bed /data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells/EGAD00010000850/processed_data/mono_meth_M_20151028_formatted.bed.gz
# --cov methylomics/QTL_TOOLS/monocytes/covariates_2.txt --permute 1000 --chunk 100 100 --exclude-phenotypes methylomics/mono_phenotypes_to_exclude.txt
# --out methylomics/QTL_TOOLS/monocytes/permutations_cov2_100_100.txt'


# assemble results
for k in covariates:
    com4 = '\'cat ' + outfolder + '/permutations_cov' + str(k + 1) + '_*_' + str(
        nchunks) + '.txt | gzip -c > ' + outfolder + '/permutations_cov' + str(k + 1) + '_full.txt.gz; '
    # "'cat methylomics/QTL_TOOLS/monocytes/permutations_cov4_*_100.txt | gzip -c > methylomics/QTL_TOOLS/monocytes/permutations_cov4_full.txt.gz; "
    com5 = '/software/R/3.4.2/bin/Rscript /home/grey2/bin/QTLtools1.1/script/runFDR_cis.R ' + \
           outfolder + '/permutations_cov' + str(
        k + 1) + '_full.txt.gz 0.05 ' + outfolder + '/permutations_cov' + str(k + 1) + '\''
    # "/software/R/3.4.2/bin/Rscript /home/grey2/bin/QTLtools1.1/script/runFDR_cis.R methylomics/QTL_TOOLS/monocytes/permutations_cov4_full.txt.gz 0.05 methylomics/QTL_TOOLS/monocytes/permutations_cov4'"
    bsubcom2 = 'bsub -o ' + outfolder + '/OUT/cov_' + str(k + 1) + '_assemble.out -e ' + outfolder + '/ERR/cov_' + str(
        k + 1) + '_assemble.err -J cov_' + str(k + 1) + '_assemble -q priority -w \'ended(\'' + phenotag + '*\')\' '
    com6 = bsubcom2 + com4 + com5
    print(com6)
    os.system(com6)
