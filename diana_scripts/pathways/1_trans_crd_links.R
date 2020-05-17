# get links of CRDs

inFolder='/Users/AvalosDiana/GitHub/IMMUNE_CELLS/diana_scripts/pathways/transCRD_associations_significant'
outFolder='/Users/AvalosDiana/GitHub/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_crds2'

trans70 = read.table(file.path(inFolder,"test.significant1FDR_trans_70.txt"),header = TRUE,stringsAsFactors=F)
trans72 = read.table(file.path(inFolder,"test.significant1FDR_trans_72.txt"),header = TRUE,stringsAsFactors=F)
trans73 = read.table(file.path(inFolder,"test.significant1FDR_trans_73.txt"),header = TRUE,stringsAsFactors=F)

translist=list(trans70, trans72, trans73)
cell_types = list('70','72','73')

for ( i in list(1,2,3) ){
  
  cell_type=cell_types[i]
  trans=translist[i][[1]]
  
  print(length(trans[,c(1)]))
  print(cell_type)
  
  trans_filt = trans[order(trans$pval),]
  trans_filt_1000=trans_filt[,c(5,10)]
  outFile=file.path(outFolder,paste0("LINKS_",cell_type,".txt"))
  write.table(file=outFile,trans_filt_1000,sep=' ')
}
