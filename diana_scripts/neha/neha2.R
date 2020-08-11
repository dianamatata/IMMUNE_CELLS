library(data.table) #useful for fread
library("readxl")
library (plyr)
library(GenomicRanges) # library(IRanges) 
library(hash) #install.packages('/Users/dianaavalos/Downloads/hash_2.2.6.1.tgz',repos = NULL)


files='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha'


################################################### characterization of enhancers

enhancers_prox <- read_excel(file.path(files,'enhancers_neha/Top50_differentially_methylated_proximal_enhancers.xlsx'))
enhancers_dist <- read_excel(file.path(files,'enhancers_neha/Top50_differentially_methylated_distal_enhancers.xlsx'))
# creation of length column
enhancers_prox$Length=enhancers_prox$End-enhancers_prox$Start
enhancers_dist$Length=enhancers_dist$End-enhancers_dist$Start
# hist plots
hist(enhancers_prox$Length, breaks=10, ylim=c(0,15), xlim=c(0,5000), col="#A42820", main ="Distribution of proximal enhancer's length", xlab='length', ylab='count')
hist(enhancers_dist$Length, breaks=10, ylim=c(0,15), xlim=c(0,2000), col="#3F5151", main ="Distribution of distal enhancer's length", xlab='length', ylab='count')
# boxplots
boxplot(enhancers_prox$Length, horizontal=FALSE, ylim=c(0,5000),frame.plot = FALSE, main ="Distribution of proximal enhancer's length", ylab='length', col="#A42820")
boxplot(enhancers_dist$Length, horizontal=FALSE, ylim=c(0,2000),frame.plot = FALSE, main ="Distribution of distal enhancer's length", ylab='length',col="#3F5151")
#adding density to hist
hist(enhancers_prox$Length, breaks=10, xlim=c(0,5000), freq = FALSE, col="#A42820", main ="Distribution of proximal enhancer's length", xlab='length', ylab='count')
lines(density(enhancers_prox$Length), col = "black",  lwd = 2)

hist(enhancers_dist$Length, breaks=10, xlim=c(0,2000), freq = FALSE, col="#3F5151", main ="Distribution of proximal enhancer's length", xlab='length', ylab='count')
lines(density(enhancers_dist$Length), col = "black",  lwd = 2)


################################################### analyse the CRDs linked

files='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha'
overlap=read.table(file.path(files,'CRD_enhancer_intersect/output_intersect_dist.txt') ,header = FALSE,stringsAsFactors=F)
colnames(overlap) = c("enhancer_chr","enhancer_start","enhancer_end","enhancer_ID","cell",
                      "CRD_chr","CRD_start","CRD_end","CRD_ID")
CRD_Length=overlap$CRD_end-overlap$CRD_start
hist(CRD_Length, breaks=10, freq = FALSE, ylim=c(0,2.3e-06), xlim=c(0,7000000), col="#AA9486", main ="Distribution of CRD's length", xlab='length', ylab='count')
lines(density(CRD_Length), col = "black",  lwd = 2)

hist(CRD_Length, ylim=c(0,60), xlim=c(0,7000000) ,breaks=10 , freq = TRUE, col="#AA9486",  main ="Distribution of CRD's length")


################################################### analysis of overlap

files='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha'
# for dist first
overlap=read.table(file.path(files,'output/output_intersect_dist.txt') ,header = FALSE,stringsAsFactors=F)
colnames(overlap) = c("enhancer_chr","enhancer_start","enhancer_end","enhancer_ID","cell",
                                   "CRD_chr","CRD_start","CRD_end","CRD_ID")
overlap$LengthEnhancer=overlap$enhancer_end-overlap$enhancer_start
# to measure overlap look at guillaume's file
enh_seq <- GRanges(seqnames=overlap$enhancer_ID,ranges=IRanges(start=overlap$enhancer_start, end=overlap$enhancer_end))
CRD_seq <- GRanges(seqnames=overlap$enhancer_ID,ranges=IRanges(start=overlap$CRD_start, end=overlap$CRD_end))
#apply(overlap, 1, function(i) min(overlap$CRD_end[i],overlap$enhancer_end[i])-max(overlap$enhancer_start[i], overlap$CRD_start[i]))
overlap$bp_overlap=mapply(function(a,b,c,d) min(c,d)-max(a,b), overlap$enhancer_start, overlap$CRD_start, overlap$CRD_end,overlap$enhancer_end)
overlap$percent_overlap=overlap$bp_overlap/overlap$LengthEnhancer*100
#plot
hist(overlap$percent_overlap, main ="overlap", xlab='length', col="#A42820")
mean(overlap$percent_overlap)
count(overlap$percent_overlap<100)
min(overlap$percent_overlap)


################################################### functions for HiC

find_HiC_overlap <- function(elem1,elem2) {# elem1 = genebed elem2 = CRDbed
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

compute_HiC_value <- function(validated,PCHiC,cell_type) {
  mapdata_validated = rep(0,nrow(mapdata))
  for(i in 1:nrow(validated)){
    currenthic = validated$queryHits[i]
    currentmap = validated$subjectHits[i]
    if (cell_type=='70'){
      currentHiCScore = PCHiC$Neu[currenthic] 
    }else if (cell_type=='72'){
      currentHiCScore = PCHiC$Mon[currenthic]
    }else if (cell_type=='73'){
      currentHiCScore = PCHiC$nCD4[currenthic] ##### depend on cell type so maybe redo plots with HiC?
    }
    if(currentHiCScore>mapdata_validated[currentmap]){
      mapdata_validated[currentmap] = currentHiCScore
    }
  }
  return(mapdata_validated)
}


###################################################  import PCHiC data
PCHiC = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks/PCHiC_peak_matrix_cutoff5.tsv')
colnames(PCHiC)[1] = "baitChr"
baitbed <- GRanges(seqnames=PCHiC$baitChr,ranges=IRanges(start=PCHiC$baitStart, end=PCHiC$baitEnd))
oebed <- GRanges(seqnames=PCHiC$oeChr,ranges=IRanges(start=PCHiC$oeStart, end=PCHiC$oeEnd))


################################################### load data

ID_enhancers_prox <- read_excel(file.path(files,'enhancers_neha/Top50_differentially_methylated_proximal_enhancers.xlsx'))[,c(1)]
ID_enhancers_dist <- read_excel(file.path(files,'enhancers_neha/Top50_differentially_methylated_distal_enhancers.xlsx'))[,c(1)]

# TRANS CRD CLUSTERS, CRD-CRD significant correlation, computed from Hi-C data, just got replaced
trans70 = read.table(file.path(files,'trans_CRD_associations_new/test.significant1FDR_trans_70.txt'),header = TRUE,stringsAsFactors=F)
trans72 = read.table(file.path(files,'trans_CRD_associations_new/test.significant1FDR_trans_72.txt'),header = TRUE,stringsAsFactors=F)
trans73 = read.table(file.path(files,'trans_CRD_associations_new/test.significant1FDR_trans_73.txt'),header = TRUE,stringsAsFactors=F)
transCRDs <- list('70' = trans70, '72' = trans70, '73' = trans73)

#CRD-gene correlations
CRDgene70=read.table(file.path(files,'CRD_gene_files/gene_CRD_mean_permutations_70.significant.txt'),stringsAsFactors=F)
CRDgene72=read.table(file.path(files,'CRD_gene_files/gene_CRD_mean_permutations_72.significant.txt'),stringsAsFactors=F)
CRDgene73=read.table(file.path(files,'CRD_gene_files/gene_CRD_mean_permutations_73.significant.txt'),stringsAsFactors=F)

CRDgene70=read.table(file.path(files,'CRD_gene_files_alt/mapping_gene_CRD_mean_70.txt'),stringsAsFactors=F)
CRDgene72=read.table(file.path(files,'CRD_gene_files_alt/mapping_gene_CRD_mean_72.txt'),stringsAsFactors=F)
CRDgene73=read.table(file.path(files,'CRD_gene_files_alt/mapping_gene_CRD_mean_73.txt'),stringsAsFactors=F)

CRDgeneList <- list('70' = CRDgene70, '72' = CRDgene72, '73' = CRDgene73)
colnames(CRDgene70)=colnames(CRDgene72)=colnames(CRDgene73) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                      "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","nominal_pval","slope","top_variant")

RES<-c("enhancer_ID","cell","CRD_chr","CRD_start","CRD_end","CRD_ID","CRD_corr_chr","CRD_corr_start","CRD_corr_end","CRD_corr_ID","corr","pval", "qval","gene_ID","pval_gene","HiC")
cell_types = list('70','72','73')



################################################### compute hic

mapdata_hic <- hash() 

for (cell_type in cell_types){
  print(cell_type)
  mapdata=data.frame(CRDgeneList[toString(cell_type)])
  mapdat=data.frame(mapdata)
  colnames(mapdat) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                       "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","nominal_pval","slope","top_variant")
  
  genebed <- GRanges(seqnames=mapdat$phenotype_ID_chr,ranges=IRanges(start=mapdat$phenotype_ID_start, end=mapdat$phenotype_ID_end))
  CRDbed <- GRanges(seqnames=mapdat$CRD_ID_chr,ranges=IRanges(start=mapdat$CRD_ID_start, end=mapdat$CRD_ID_end))
  validated_gene_CRD = find_HiC_overlap(genebed,CRDbed)
  mapdat$HIC=compute_HiC_value(validated_gene_CRD,PCHiC,cell_type)
  mapdata_hic[[toString(cell_type)]] <- mapdat
}

###
cal
# TO CHANGE 
ID_enhancers=ID_enhancers_prox
CRD_overlapping_list=read.table(file.path(files,'CRD_enhancer_intersect/output_intersect_prox.txt') ,header = FALSE,stringsAsFactors=F)
outputfile='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha/results_prox_hic'

ID_enhancers=ID_enhancers_dist
CRD_overlapping_list=read.table(file.path(files,'CRD_enhancer_intersect/output_intersect_dist.txt') ,header = FALSE,stringsAsFactors=F)
outputfile='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha/results_dist1_alt'


###################################################  compute
colnames(CRD_overlapping_list) = c("enhancer_chr","enhancer_start","enhancer_end","enhancer_ID","cell",
                      "CRD_chr","CRD_start","CRD_end","CRD_ID")

for (idx in 1:lengths(ID_enhancers)){
  ID_enhancer = ID_enhancers$ID[idx]
  CRD_overlapping_enhancers_list=dplyr::filter(CRD_overlapping_list, grepl(ID_enhancer, enhancer_ID))
  CRD_overlap=CRD_overlapping_enhancers_list$CRD_ID
  #print(paste(idx,'CRD_overlap',CRD_overlap,'ID_enhancer ', ID_enhancer))
  # for each CRD overlapping in each cell type
  if (length(CRD_overlap)>0){
    for (idx_crd in 1:length(CRD_overlap)){
      CRDo=CRD_overlap[idx_crd] # CRD overlapping i
      cell_type=CRD_overlapping_enhancers_list$cell[idx_crd]
      trans=transCRDs[[toString(cell_type)]]  
      CRDgenes=mapdata_hic[[toString(cell_type)]]
      #print(paste(idx_crd,CRDo,cell_type)) 
      
      a=dplyr::filter(trans, grepl(CRDo, id1)) [,c(7,8,9,10,11,12,13)]
      b=dplyr::filter(trans, grepl(CRDo, id2)) [,c(2,3,4,5,11,12,13)]
      colnames(b)=colnames(a)
      crds_correlated=list(a,b)
      crds_correlated_ids=append(unique(b$id1),unique(a$id2))
      crds_correlated_df <- ldply (crds_correlated, data.frame) #list to df
      
      if (length(crds_correlated_ids)>0){
        for (idx_corr in 1:length(crds_correlated_ids)){
          genes=dplyr::filter(CRDgenes, grepl(crds_correlated_ids[idx_corr], CRD_ID))
          if(nrow(genes)>0){
            # print(paste(" idx_corr ",idx_corr,crds_correlated_ids[idx_corr]," ngenes ", toString(nrow(genes))) )
            for (idx_gene in 1:nrow(genes)){ # if many genes
              #print(genes$V8[idx_gene])#crd linked to gene
              corr_crds=dplyr::filter(crds_correlated_df, grepl(crds_correlated_ids[idx_corr], id2))
              row=unlist(list(ID_enhancer,(dplyr::filter(CRD_overlapping_list, grepl(CRDo, CRD_ID)) [1,c(5,6,7,8,9)]),corr_crds, genes[idx_gene,c(1,12,15)]))
              #print(row)
              ROW<- ldply (row, data.frame) 
              c=vapply(ldply (row, data.frame)[,c(2)], paste, collapse = ", ", character(1L))
              row_df<- t(ldply (c, data.frame))
              colnames(row_df)=c("enhancer_ID","cell","CRD_chr","CRD_start","CRD_end","CRD_ID","CRD_corr_chr","CRD_corr_start",
                                 "CRD_corr_end","CRD_corr_ID","corr","pval", "qval","gene_ID","nominal_pval_gene","hic")
              RES=rbind(RES,row_df)
            }
          }
        }
      }
    }
  }
}
RESULT=as.data.frame(RES)
table(RESULT$cell)

# CHANGE
RES_prox=as.data.frame(RES)

RES_dist=as.data.frame(RES)

#

RES_prox1=as.data.frame(RES_prox)
RES_prox1$hic <- as.numeric(as.character(RES_prox1$hic))
RES_prox_hic_cutoff=filter(RES_prox1, hic > 5)
table(RES_prox_hic_cutoff$cell)
write.table(file='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha/results_prox_hic',RES_prox1,quote=F, row.names=F,sep=' ')
write.table(file='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha/results_prox_hic_cutoff5',RES_prox_hic_cutoff,quote=F, row.names=F,sep=' ')
table(RES_prox1$cell)
table(RES_prox_hic_cutoff$cell)


RES_dist1=as.data.frame(RES_dist)
RES_dist1$hic <- as.numeric(as.character(RES_dist1$hic))
RES_dist_hic_cutoff=filter(RES_dist1, hic > 5)
table(RES_prox_hic_cutoff$cell)
write.table(file='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha/results_dist_hic',RES_dist1,quote=F, row.names=F,sep=' ')
write.table(file='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha/results_dist_hic_cutoff5',RES_dist_hic_cutoff,quote=F, row.names=F,sep=' ')
table(RES_dist1$cell)
table(RES_dist_hic_cutoff$cell)



