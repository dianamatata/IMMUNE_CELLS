# use the libraries available to see the GO of our list of genes
# Clean environment
rm(list=ls())
gc()library(gprofiler2)

sources_in=c("GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","TF","MIRNA","HPA","CORUM","HP")
sources_in=c("GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","HPA","CORUM","HP")
path='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_genes'
path2='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_pathways'
path='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:TRHs/trans_hubs_genes'
path2='/Users/dianaavalos/Programming/A_CRD_plots/trans_files/7_CRD_Trans:TRHs/trans_hubs_pathways'

files <- Sys.glob(file.path(path, "*genesonly*"))

# f='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_genes/genesonly_73_cluster2.csv'
for(f in files){
  genelist=list(read.table(f, header = FALSE, sep = "", dec = "."))
  gostres <- gost(query = genelist, organism = "hsapiens", sources=sources_in, user_threshold=1e-5)
  # head(gostres$result)
  # p <- gostplot(gostres, capped = FALSE, interactive = TRUE)
  # p
  if (!is.null(gostres$result)){
    print(unlist(strsplit(f, "genesonly_"))[2])
    print(length(genelist))
    gostres$result=gostres$result[order(gostres$result$p_value),] # order by pval
    outfile=file.path(path2, unlist(strsplit(f, "genesonly_"))[2])
    fwrite(gostres$result, file =outfile, sep='\t')
    # if(length(gostres$result$term_id)){
    #   outfile2=paste0(path2, '/G_',unlist(strsplit(f, "genesonly_"))[2])
    #   fwrite(list(gostres$result$term_id), file =outfile2, sep='\t')
    # }
  }
}

#### second method

library("GOstats")
library(GOstats)


## third
