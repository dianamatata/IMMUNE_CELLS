library(tidyverse)
library(data.table)

files <- Sys.glob(file.path("for_RNA2", "*gene.rpkm.bed.gz"))
for (filename in files){
  mydat = fread(filename)
  mydat.quant = mydat[,-c(1:6)]
  goodlines = rowSums(mydat.quant == 0) <= 20 # remove genes with low coverage. length(goodlines[goodlines == TRUE])
  # remove paired-end samples. to find them: cd ~/scratch/Blueprint/RNA_seq/EGAD00001002675% ls | grep 'paired' | cut -c1-8 | uniq
  samples_to_remove = c("S001C2B4","S001T511","S001YW11","S0020M11","S003JHB5","S006UKB2","S00DKCB3","S00JV3B4","S00JYYB2")
  colnames(mydat)[-c(1:6)] = unlist(lapply(colnames(mydat)[-c(1:6)],substring,3,10))
  mydat.filtered = mydat[goodlines ,!colnames(mydat) %in% samples_to_remove, with=FALSE]
  colnames(mydat.filtered)[-c(1:6)] = unlist(lapply(colnames(mydat.filtered)[-c(1:6)],substring,1,6))
  output_file=str_replace(filename, 'quantification_', 'quantification_filtered_')
  write.table(file='qtltools_quantification_filtered75.gene.rpkm.bed',mydat.filtered,quote=F,row.names=F,sep='\t')
}

