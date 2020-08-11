# library(qvalue)
# library(ggplot2)
install.packages("gprofiler2")
library(data.table) #useful for fread
library(gdata) #for reading excel files
library(GenomicRanges) # for GRanges function in HIC data, findOverlaps defined from package "IRanges"
# library(tidyverse)

### CRDs doverlapping differentially methylated enhancers
outputs='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha/output'

dist70_CRDs=read.table(file.path(outputs,'output_intersect_dist_70.txt'),header = FALSE,stringsAsFactors=F)[,c(8)]
dist72_CRDs=read.table(file.path(outputs,'output_intersect_dist_72.txt'),header = FALSE,stringsAsFactors=F)[,c(8)]
dist73_CRDs=read.table(file.path(outputs,'output_intersect_dist_73.txt'),header = FALSE,stringsAsFactors=F)[,c(8)]
prox70_CRDs=read.table(file.path(outputs,'output_intersect_prox_70.txt'),header = FALSE,stringsAsFactors=F)[,c(8)]
prox72_CRDs=read.table(file.path(outputs,'output_intersect_prox_72.txt'),header = FALSE,stringsAsFactors=F)[,c(8)]
prox73_CRDs=read.table(file.path(outputs,'output_intersect_prox_73.txt'),header = FALSE,stringsAsFactors=F)[,c(8)]

# TRANS CRD CLUSTERS, CRD-CRD significant correlation
path='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/transCRD_associations_significant'
trans70 = read.table(file.path(path,'test.significant1FDR_trans_70.txt'),header = TRUE,stringsAsFactors=F)
trans72 = read.table(file.path(path,'test.significant1FDR_trans_72.txt'),header = TRUE,stringsAsFactors=F)
trans73 = read.table(file.path(path,'test.significant1FDR_trans_73.txt'),header = TRUE,stringsAsFactors=F)

translist=list(trans70, trans72, trans73)
cell_types = list('70','72','73')
CRDs_dist  = list(dist70_CRDs, dist72_CRDs, dist73_CRDs)
CRDs_prox  = list(prox70_CRDs, prox72_CRDs, prox73_CRDs)


### filter CRDs, find the correlations in the trans files

# filter for id1 and id2 columns
crds_filtered=with(trans70, trans70[ grepl(paste(dist70_CRDs, collapse = "|"), id1) | grepl(paste(dist70_CRDs, collapse = "|"), id2), ])
crds_list=append(unique(crds_filtered$id1),unique(crds_filtered$id2))
length(crds_list)
#look at cis genes associated
CRD_genes_list_70 = read.table("/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/cisCRD-gene_associations/mapping_gene_CRD_mean_ALL_70.txt", header=FALSE, sep=' ')
dplyr::filter(CRD_genes_list_70, grepl(paste(crds_list, collapse = "|"), V8))



# look if these CRDs are connected 


for ( i in list(1,2,3) ){
  cell_type=cell_types[i]
  trans=translist[i][[1]]
  
  
}

  
  
  
  
  


# EGAD00001002673_CRDs_info.MOD1.NRE2.txt.gz for CRD name chr start end



# concatenate the 3 CRDs info
CRDs_all_cells = data.frame()
for(i in list(70, 72,73)){
  CRDinfo = fread(file.path(path, sprintf('EGAD000010026%d_CRDs_info.MOD1.NRE2.txt.gz',i)),header=F)
  CRDinfo$cell=i
  CRDs_all_cells=rbind(CRDs_all_cells,CRDinfo)
  colnames(CRDs_all_cells) = c("CRD_ID","CRD_chr","CRD_start","CRD_end","Cell")
}


## import files
# annotations
path='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha'
protein_coding_genes=as.data.frame(rtracklayer::import(file.path(path, 'gencode.v19.chr_patch_hapl_scaff.annotation.gtf')))
long_nc_RNA_genes <- rtracklayer::import(file.path(path, 'gencode.v19.long_noncoding_RNAs.gtf'))
long_nc_RNA_genes=as.data.frame(long_nc_RNA_genes)
genelist = c(protein_coding_genes,long_nc_RNA_genes)

# CRDinfo
CRDinfo = fread(file.path(path, 'EGAD00001002673_CRDs_info.MOD1.NRE2.txt.gz'),header=F)
a = rep(0., length(CRDinfo$V2))
for(i in 1:length(CRDinfo$V2)) {a[i]=paste0("chr", toString(CRDinfo$V2[i])) }
CRDinfo$V2=a
CRDinfo <- CRDinfo[, c(2, 3, 4, 1)]
write.table(file=file.path(path, 'CRDinfo73.txt'),CRDinfo,quote=F,row.names=F,sep='\t', col.names=FALSE)

# Neha
meth_dist = read.xls (file.path(path, 'Top50_differentially_methylated_distal_enhancers.xlsx'), sheet = 1, header = TRUE)
meth_prox = read.xls (file.path(path, 'Top50_differentially_methylated_proximal_enhancers.xlsx'), sheet = 1, header = TRUE)

meth_dist <- meth_dist[, c(2, 3, 4, 1)]
meth_prox <- meth_prox[, c(2, 3, 4, 1)]
write.table(file=file.path(path, 'meth_dist.txt'),meth_dist,quote=F,row.names=F,sep='\t', col.names=FALSE)
write.table(file=file.path(path, 'meth_prox.txt'),meth_prox,quote=F,row.names=F,sep='\t', col.names=FALSE)

########

# cat CRDinfo73.txt | sort -V -k1,1 -k2,2n > CRDinfo73s.txt

mapdata = read.table('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/mapping_gene_CRD_mean_ALL_70.txt',stringsAsFactors=F)
colnames(mapdata) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                      "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","nominal_pval","slope","top_variant")

# remove chr7 to 7
meth_dist$Chromosome=by(meth_dist, 1:nrow(meth_dist), function(row) unlist(strsplit(as.character(row$Chromosome),"r"))[2] )
meth_prox$Chromosome=by(meth_prox, 1:nrow(meth_prox), function(row) unlist(strsplit(as.character(row$Chromosome),"r"))[2] )


# HiC
PCHiC = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/PCHiC_peak_matrix_cutoff5.tsv')
colnames(PCHiC)[1] = "baitChr"
baitbed <- GRanges(seqnames=PCHiC$baitChr,ranges=IRanges(start=PCHiC$baitStart, end=PCHiC$baitEnd))
oebed <- GRanges(seqnames=PCHiC$oeChr,ranges=IRanges(start=PCHiC$oeStart, end=PCHiC$oeEnd))
genebed <- GRanges(seqnames=mapdata$phenotype_ID_chr,ranges=IRanges(start=mapdata$phenotype_ID_start, end=mapdata$phenotype_ID_end))
CRDbed <- GRanges(seqnames=mapdata$CRD_ID_chr,ranges=IRanges(start=mapdata$CRD_ID_start, end=mapdata$CRD_ID_end))

methdistbed <- GRanges(seqnames=as.character(meth_dist$Chromosome),ranges=IRanges(start=meth_dist$Start, end=meth_dist$End))
methproxbed <- GRanges(seqnames=as.character(meth_prox$Chromosome),ranges=IRanges(start=meth_prox$Start, end=meth_prox$End))

validated = find_HiC_overlap(genebed,CRDbed)
validated = find_HiC_overlap(genebed,methdistbed)
validated = find_HiC_overlap(genebed,methproxbed)

hic_validated = rep(1,nrow(PCHiC))

for(i in 1:nrow(validated)){
  currenthic = validated$queryHits[i]
  currentpval = mapdata$nominal_pval[validated$subjectHits[i]]
  if(currentpval<hic_validated[currenthic]){
    hic_validated[currenthic] = currentpval
  }
}

mapdata_validated = rep(0,nrow(mapdata))

for(i in 1:nrow(validated)){
  currenthic = validated$queryHits[i]
  currentmap = validated$subjectHits[i]
  currentHiCScore = PCHiC$nCD4[currenthic]
  if(currentHiCScore>mapdata_validated[currentmap]){
    mapdata_validated[currentmap] = currentHiCScore
  }
}








