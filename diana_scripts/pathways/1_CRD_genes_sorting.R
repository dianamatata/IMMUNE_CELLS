# GENES TO CRD

# SCRIPT OBSOLETE. we get the TRHs from /Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/trans/TRHIndex_3_cell_types.R

# For each CRD, look at the genes they are regulating in cis. (careful, how do we know if it is the gene regulating the CRD or the CRD regulating the gene, which direction?)

mapdataNEU = read.table('/Users/dianaavalos/Programming/pathways/mapping_gene_CRD_mean_ALL_70.txt',stringsAsFactors=F)
mapdataMON = read.table('/Users/dianaavalos/Programming/pathways/mapping_gene_CRD_mean_ALL_72.txt',stringsAsFactors=F)
mapdataTCL = read.table('/Users/dianaavalos/Programming/pathways/mapping_gene_CRD_mean_ALL_73.txt',stringsAsFactors=F)
mapdatas=list(mapdataNEU, mapdataMON, mapdataTCL)
cell_types = list('NEU','MON','TCL')

for ( i in list(1,2,3) ){
  
  cell_type=cell_types[i]
  mapdata=mapdatas[i][[1]]
  
  colnames(mapdata) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                        "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","nominal_pval","slope","top_variant")
  
  # next step: look up all the genes associated with one CRD, put them in a list
  CRD_unique_IDs=unique(mapdata$CRD_ID) #CRD unique list
  out_file=paste0("/Users/dianaavalos/Programming/pathways/CRD-genes-associated_", cell_type, ".txt") # filename
  cat('',file=out_file,sep=" ") # create file
  for (CRD_i in CRD_unique_IDs){
    associated_genes=mapdata[which(mapdata$CRD_ID == CRD_i),c(1:5)]
    if (length(associated_genes$phenotype_ID)>1){
      cat(paste0(CRD_i, ' '),file=out_file,append=TRUE)
      cat(print(mapdata[which(mapdata$CRD_ID == CRD_i),c(1)]),file=out_file,append=TRUE)
      cat("\n",file=out_file,append=TRUE)

    }
  }
}


# TRANS CRD CLUSTERS
path='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/transCRD_associations_significant'
trans70 = read.table(file.path(path,'test.significant1FDR_trans_70.txt'),header = TRUE,stringsAsFactors=F)
trans72 = read.table(file.path(path,'test.significant1FDR_trans_72.txt'),header = TRUE,stringsAsFactors=F)
trans73 = read.table(file.path(path,'test.significant1FDR_trans_73.txt'),header = TRUE,stringsAsFactors=F)

translist=list(trans70, trans72, trans73)

for ( i in list(1,2,3) ){
  
  cell_type=cell_types[i]
  trans=translist[i][[1]]
  
  print(length(trans[,c(1)]))
  print(cell_type)
  
  trans_filt = trans[order(trans$pval),]
#trans_filt = trans_filt[1:1000, ]
trans_filt_1000=trans_filt[,c(5,10)]
write.table(file='/Users/dianaavalos/Programming/pathways/LINKS_70.txt',trans_filt_1000,quote=F, row.names=F,sep=' ')

# export and use python code to have all the clusters in some file
# run /Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/trans/links_hubs.py to get the CRDs in the same clusters
# then for each cell type, take all the genes included in the same cluster of CRDs and list thm in a text file taking mapdata[which(mapdata$CRD_ID == CRD_i),c(1:5)]

