cat gene_CRD_mean_chunk*.txt | gzip -c > gene_CRD_mean_permutations_full.txt.gz
/software/R/3.4.2/bin/Rscript /home/grey2/bin/QTLtools1.1/script/runFDR_cis.R gene_CRD_mean_permutations_full.txt.gz 0.01 gene_CRD_mean_permutations
