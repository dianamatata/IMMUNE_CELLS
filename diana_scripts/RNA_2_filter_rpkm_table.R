library(tidyverse)
library(data.table)
mydat = fread('qtltools_quantification.gene.rpkm.bed.gz')
mydat.quant = mydat[,-c(1:6)]
goodlines = rowSums(mydat.quant == 0) <= 20
#remove paired-end samples
samples_to_remove = c("S001C2B4","S001T511","S001YW11","S0020M11","S003JHB5","S006UKB2","S00DKCB3","S00JV3B4","S00JYYB2")
colnames(mydat)[-c(1:6)] = unlist(lapply(colnames(mydat)[-c(1:6)],substring,3,10))
mydat.filtered = mydat[goodlines ,!colnames(mydat) %in% samples_to_remove, with=FALSE]
colnames(mydat.filtered)[-c(1:6)] = unlist(lapply(colnames(mydat.filtered)[-c(1:6)],substring,1,6))
write.table(file='qtltools_quantification_filtered.gene.rpkm.bed',mydat.filtered,quote=F,row.names=F,sep='\t')

