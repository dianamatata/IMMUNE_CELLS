library(data.table) #useful for fread
library("readxl")
library (plyr)
library(GenomicRanges) # for GRanges function in HIC data, findOverlaps defined from package "IRanges"


find_HiC_overlap <- function(elem1,elem2) {
  # elem1 = genebed elem2 = CRDbed
  
  #fwd
  x = findOverlaps(baitbed,elem1)
  y = findOverlaps(oebed,elem2)
  tmp = rbind(as.data.frame(x),as.data.frame(y))  # rbind() function combines vector, matrix or data frame by rows.
  validated.fwd = tmp[which(duplicated(tmp)),]
  
  #bwd
  x = findOverlaps(baitbed,elem2)
  y = findOverlaps(oebed,elem1)
  tmp = rbind(as.data.frame(x),as.data.frame(y))
  validated.bwd = tmp[which(duplicated(tmp)),]
  
  validated = unique(rbind(validated.fwd,validated.bwd))
  return(validated)
}

files='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha'
ID_enhancers_prox <- read_excel(file.path(files,'Top50_differentially_methylated_proximal_enhancers.xlsx'))[,c(1)]
ID_enhancers_dist <- read_excel(file.path(files,'Top50_differentially_methylated_distal_enhancers.xlsx'))[,c(1)]

# TRANS CRD CLUSTERS, CRD-CRD significant correlation
path='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/transCRD_associations_significant'
trans70 = read.table(file.path(path,'test.significant1FDR_trans_70.txt'),header = TRUE,stringsAsFactors=F)
trans72 = read.table(file.path(path,'test.significant1FDR_trans_72.txt'),header = TRUE,stringsAsFactors=F)
trans73 = read.table(file.path(path,'test.significant1FDR_trans_73.txt'),header = TRUE,stringsAsFactors=F)
transCRDs <- list('70' = trans70, '72' = trans70, '73' = trans73)

#CRD-gene correlations
path2='/Users/dianaavalos/Programming/Hi-C_correlated_peaks'
CRDgene70=read.table(file.path(path2,'mapping_gene_CRD_mean_ALL_70.txt'),stringsAsFactors=F)
CRDgene72=read.table(file.path(path2,'mapping_gene_CRD_mean_ALL_72.txt'),stringsAsFactors=F)
CRDgene73=read.table(file.path(path2,'mapping_gene_CRD_mean_ALL_73.txt'),stringsAsFactors=F)
CRDgeneList <- list('70' = CRDgene70, '72' = CRDgene72, '73' = CRDgene73)
colnames(CRDgene70)=colnames(CRDgene72)=colnames(CRDgene73) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                                                                "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","nominal_pval","slope","top_variant")

cell_types = list('70','72','73')

# TO CHANGE
ID_enhancers=ID_enhancers_prox
CRD_overlapping_list=read.table(file.path(files,'output/output_intersect_prox.txt') ,header = FALSE,stringsAsFactors=F)
outputfile='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha/results_prox'

ID_enhancers=ID_enhancers_dist
CRD_overlapping_list=read.table(file.path(files,'output/output_intersect_dist.txt') ,header = FALSE,stringsAsFactors=F)
outputfile='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha/results_dist'


colnames(CRD_overlapping_list) = c("enhancer_chr","enhancer_start","enhancer_end","enhancer_ID","cell",
                                   "CRD_chr","CRD_start","CRD_end","CRD_ID")

for (idx in 1:lengths(ID_enhancers)){
  ID_enhancer = ID_enhancers$ID[idx]
  CRD_overlapping_enhancers_list=dplyr::filter(CRD_overlapping_list, grepl(ID_enhancer, enhancer_ID))
  CRD_overlap=CRD_overlapping_enhancers_list$CRD_ID
  # for each CRD overlapping in each cell type
  if (length(CRD_overlap)>0){
    for (idx_crd in 1:length(CRD_overlap)){
      CRDo=CRD_overlap[idx_crd] # CRD overlapping i
      cell_type=CRD_overlapping_enhancers_list$cell[idx_crd]
      
    }
  }
}

mapdata = CRDgene70

PCHiC = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/PCHiC_peak_matrix_cutoff5.tsv')
colnames(PCHiC)[1] = "baitChr"
CRD_overlapping_list
baitbed <- GRanges(seqnames=PCHiC$baitChr,ranges=IRanges(start=PCHiC$baitStart, end=PCHiC$baitEnd))
oebed <- GRanges(seqnames=PCHiC$oeChr,ranges=IRanges(start=PCHiC$oeStart, end=PCHiC$oeEnd))
genebed <- GRanges(seqnames=mapdata$phenotype_ID_chr,ranges=IRanges(start=mapdata$phenotype_ID_start, end=mapdata$phenotype_ID_end))
CRDbed <- GRanges(seqnames=mapdata$CRD_ID_chr,ranges=IRanges(start=mapdata$CRD_ID_start, end=mapdata$CRD_ID_end))

methdistbed <- GRanges(seqnames=as.character(meth_dist$Chromosome),ranges=IRanges(start=meth_dist$Start, end=meth_dist$End))
methproxbed <- GRanges(seqnames=as.character(meth_prox$Chromosome),ranges=IRanges(start=meth_prox$Start, end=meth_prox$End))

validated = find_HiC_overlap(genebed,CRDbed)
validated = find_HiC_overlap(genebed,CRD_overlapping_list)

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


