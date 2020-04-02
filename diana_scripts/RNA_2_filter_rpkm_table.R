library(tidyverse)
library(data.table)

for (cell_type in list('EGAD00001002671', 'EGAD00001002674' , 'EGAD00001002675')){
	files <- Sys.glob(file.path("/home/users/a/avalosma/scratch/RNA_quantify_indep", paste0("*", cell_type, "*gene.rpkm.bed")))
	print (length(files))
	common_data=fread(files[1])[,c(1:6)]
	for (filename in files){
		# print(filename)
		mydat = fread(filename)
		colnames(mydat)[7] = unlist(strsplit( unlist(strsplit( unlist(strsplit(colnames(mydat)[-c(1:6)], "[/]"))[10] , "[_]"))[1],  "[.]"))[1] # rename sample
		common_data$name=mydat[,c(7)] # add it to the common data table of the cell type
		names(common_data)[names(common_data) == "name"] <- colnames(mydat)[7] # change the name of the column with the right sample
	}

	common_data.quant = common_data[,-c(1:6)]
	goodlines = rowSums(common_data.quant == 0) <= 20 # remove genes with low coverage. length(goodlines[goodlines == TRUE])
	# remove paired-end samples. to find them: cd ~/scratch/Blueprint/RNA_seq/EGAD00001002675% ls | grep 'paired' | cut -c1-8 | uniq
	samples_to_remove = c("S001C2B4","S001T511","S001YW11","S0020M11","S003JHB5","S006UKB2","S00DKCB3","S00JV3B4","S00JYYB2")
	common_data.filtered = mydat[goodlines ,!colnames(common_data) %in% samples_to_remove, with=FALSE]
	output_file = paste0("/home/users/a/avalosma/scratch/RNA_rpkm/", cell_type, "_quantification_filtered.gene.rpkm.bed") 
	write.table(file=output_file,common_data.filtered,quote=F,row.names=F,sep='\t')
}
