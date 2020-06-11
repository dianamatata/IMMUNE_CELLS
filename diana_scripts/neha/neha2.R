library(data.table) #useful for fread
library("readxl")
library (plyr)


files='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/neha'

ID_enhancers_prox <- read_excel(file.path(files,'Top50_differentially_methylated_proximal_enhancers.xlsx'))[,c(1)]
ID_enhancers_dist <- read_excel(file.path(files,'Top50_differentially_methylated_distal_enhancers.xlsx'))[,c(1)]

# TRANS CRD CLUSTERS, CRD-CRD significant correlation, computed from Hi-C data
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

RES<-c("enhancer_ID","cell","CRD_chr","CRD_start","CRD_end","CRD_ID","CRD_corr_chr","CRD_corr_start","CRD_corr_end","CRD_corr_ID","corr","pval", "qval","gene_ID","pval_gene")
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
  # print(paste(idx,'CRD_overlap',CRD_overlap,'ID_enhancer ', ID_enhancer))
  # for each CRD overlapping in each cell type
  if (length(CRD_overlap)>0){
    for (idx_crd in 1:length(CRD_overlap)){
      CRDo=CRD_overlap[idx_crd] # CRD overlapping i
      cell_type=CRD_overlapping_enhancers_list$cell[idx_crd]
      trans=transCRDs[[toString(cell_type)]]  
      CRDgenes=CRDgeneList[[toString(cell_type)]] 
      print(paste(idx_crd,CRDo,cell_type)) 
      
      a=dplyr::filter(trans, grepl(CRDo, id1)) [,c(7,8,9,10,11,12,13)]
      b=dplyr::filter(trans, grepl(CRDo, id2)) [,c(2,3,4,5,11,12,13)]
      colnames(b)=colnames(a)
      crds_correlated=list(a,b)
      crds_correlated_ids=append(unique(b$id1),unique(a$id2))
      crds_correlated_df <- ldply (crds_correlated, data.frame) #list to df
      
      if (length(crds_correlated_ids)>0){
        for (idx_corr in 1:length(crds_correlated_ids)){
          genes=dplyr::filter(CRDgenes, grepl(crds_correlated_ids[idx_corr], V8))
          if(nrow(genes)>0){
            # print(paste(" idx_corr ",idx_corr,crds_correlated_ids[idx_corr]," ngenes ", toString(nrow(genes))) )
            for (idx_gene in 1:nrow(genes)){ # if many genes
              #print(genes$V8[idx_gene])#crd linked to gene
              corr_crds=dplyr::filter(crds_correlated_df, grepl(crds_correlated_ids[idx_corr], id2))
              row=unlist(list(ID_enhancer,(dplyr::filter(CRD_overlapping_list, grepl(CRDo, CRD_ID)) [1,c(5,6,7,8,9)]),corr_crds, genes[idx_gene,c(1,12)]))
              #print(row)
              ROW<- ldply (row, data.frame) 
              c=vapply(ldply (row, data.frame)[,c(2)], paste, collapse = ", ", character(1L))
              row_df<- t(ldply (c, data.frame))
              colnames(row_df)=c("enhancer_ID","cell","CRD_chr","CRD_start","CRD_end","CRD_ID","CRD_corr_chr","CRD_corr_start",
                                 "CRD_corr_end","CRD_corr_ID","corr","pval", "qval","gene_ID","nominal_pval_gene")
              RES=rbind(RES,row_df)
            }
          }
        }
      }
    }
  }
}

write.table(file=outputfile,RES,quote=F, row.names=F,sep=' ')


