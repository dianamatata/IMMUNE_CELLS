#!/bin/bash

K=100
BIN=/data/unige/funpopgen/odelanea/SGX/V2/binaries/qtltools/bin/QTLtools

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002675_v3.0/EGAD00001002675_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/70_vs_72_ALL.ALLchr.mean.txt.gz

for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/70_vs_72/gene_CRD_mean_chunk$k
        #echo "$BIN cis --vcf $MOD --bed $RNA --permute 200 --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002671_v3.0/EGAD00001002671_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/73_vs_72_ALL.ALLchr.mean.txt.gz

for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/73_vs_72/gene_CRD_mean_chunk$k
        #echo "$BIN cis --vcf $MOD --bed $RNA --permute 200 --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002675_v3.0/EGAD00001002675_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/70_vs_73_ALL.ALLchr.mean.txt.gz

for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/70_vs_73/gene_CRD_mean_chunk$k
        #echo "$BIN cis --vcf $MOD --bed $RNA --permute 200 --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002674_v3.0/EGAD00001002674_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/72_vs_73_ALL.ALLchr.mean.txt.gz

for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/72_vs_73/gene_CRD_mean_chunk$k
        #echo "$BIN cis --vcf $MOD --bed $RNA --permute 200 --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002674_v3.0/EGAD00001002674_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/72_vs_70_ALL.ALLchr.mean.txt.gz

for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/72_vs_70/gene_CRD_mean_chunk$k
        #echo "$BIN cis --vcf $MOD --bed $RNA --permute 200 --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002671_v3.0/EGAD00001002671_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/73_vs_70_ALL.ALLchr.mean.txt.gz

for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/73_vs_70/gene_CRD_mean_chunk$k
        #echo "$BIN cis --vcf $MOD --bed $RNA --permute 200 --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done

#exit

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002675_v3.0/EGAD00001002675_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/70_vs_72_ALL.ALLchr.mean.txt.gz
TH=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/70_vs_72/gene_CRD_mean_permutations.thresholds.txt
for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/70_vs_72/mapping_gene_CRD_mean_chunk$k
        #echo "$BIN cis --vcf $MOD --bed $RNA --nominal $TH --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002671_v3.0/EGAD00001002671_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/73_vs_72_ALL.ALLchr.mean.txt.gz
TH=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/73_vs_72/gene_CRD_mean_permutations.thresholds.txt
for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/73_vs_72/mapping_gene_CRD_mean_chunk$k
        echo "$BIN cis --vcf $MOD --bed $RNA --nominal $TH --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002674_v3.0/EGAD00001002674_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/72_vs_70_ALL.ALLchr.mean.txt.gz
TH=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/72_vs_70/gene_CRD_mean_permutations.thresholds.txt
for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/72_vs_70/mapping_gene_CRD_mean_chunk$k
        echo "$BIN cis --vcf $MOD --bed $RNA --nominal $TH --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002671_v3.0/EGAD00001002671_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/73_vs_70_ALL.ALLchr.mean.txt.gz
TH=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/73_vs_70/gene_CRD_mean_permutations.thresholds.txt
for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/73_vs_70/mapping_gene_CRD_mean_chunk$k
        echo "$BIN cis --vcf $MOD --bed $RNA --nominal $TH --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002675_v3.0/EGAD00001002675_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/70_vs_73_ALL.ALLchr.mean.txt.gz
TH=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/70_vs_73/gene_CRD_mean_permutations.thresholds.txt
for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/70_vs_73/mapping_gene_CRD_mean_chunk$k
        #echo "$BIN cis --vcf $MOD --bed $RNA --nominal $TH --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done

RNA=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/QTL_TOOLS/EGAD00001002674_v3.0/EGAD00001002674_RNA.PC10.bed.gz
MOD=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/72_vs_73_ALL.ALLchr.mean.txt.gz
TH=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/72_vs_73/gene_CRD_mean_permutations.thresholds.txt
for k in $(seq 1 $K); do
        OUT=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/mapping_aCRD_gene/72_vs_73/mapping_gene_CRD_mean_chunk$k
        #echo "$BIN cis --vcf $MOD --bed $RNA --nominal $TH --chunk $k $K --out $OUT\.txt" | bsub -o $OUT\.log
done
