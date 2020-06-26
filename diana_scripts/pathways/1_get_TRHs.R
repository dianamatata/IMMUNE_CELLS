
library(qvalue)
library(ggplot2)
library(gplots)
library(data.table)
library(tidyverse)
library(igraph)


# keep DATA when qvalue < 0.01
DATA_70 = as.data.frame(data.table::fread("/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/EGAD00001002670_ALL.ALLchr.txt.gz", head=FALSE, stringsAsFactors=FALSE))
DATA_72 = as.data.frame(data.table::fread("/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/EGAD00001002672_ALL.ALLchr.txt.gz", head=FALSE, stringsAsFactors=FALSE))
DATA_73 = as.data.frame(data.table::fread("/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/EGAD00001002673_ALL.ALLchr.txt.gz", head=FALSE, stringsAsFactors=FALSE))

DATA=DATA_70 ######## cell type
colnames(DATA) = c("idx1","chr1","start1","end1","id1","idx2","chr2","start2","end2","id2","corr","pval")
Q = qvalue(DATA$pval)
DATA$qval = Q$qvalues
DATAs = DATA[Q$qvalue < 0.01, ]

###### CALL TRHs using greedy algorithm
LINKS.TRH = data.frame(from=DATAs$id1, to=DATAs$id2, weigth=ifelse(DATAs$corr>0, 1, 2), stringsAsFactors=FALSE)
NODES.TRH = data.frame(id=names(table(c(LINKS.TRH$from, LINKS.TRH$to))), chr=matrix(unlist(strsplit(names(table(c(LINKS.TRH$from, LINKS.TRH$to))), split="_")), ncol=3, byrow=TRUE)[, 1], stringsAsFactors=FALSE)
DATAnet.TRH = graph_from_data_frame(d=LINKS.TRH, vertices=NODES.TRH, directed=FALSE)
# do not plot takes a long time and aborts session #plot(DATAnet.TRH, vertex.label = V(DATAnet.TRH)$name)
communities  = fastgreedy.community(DATAnet.TRH) # https://www.rdocumentation.org/packages/igraph/versions/0.4.1/topics/fastgreedy.community
# this function tries to find dense subgraph, also called communities in graphs via directly optimizing a modularity score.
N_TRH = max(communities$membership)
N_TRH_TOSHOW = 50

# number of TRH
length(unique(communities$membership))
# number of CRDs involved
length(communities$names)
COM = data.frame(communities$names, communities$membership)
COM = COM[order(communities$membership),]

path='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/trans'
######## cell type
write.table(file=file.path(path,'TRH_70.txt'),COM,quote=F,row.names=F,sep='\t')