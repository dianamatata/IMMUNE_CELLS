# Goal: plot 50 TRH index for 3 cell types

library(qvalue)
library(ggplot2)
library(gplots)
library(data.table)
library(tidyverse)


# keep DATA when qvalue < 0.01
DATA_70 = as.data.frame(data.table::fread("/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/EGAD00001002670_ALL.ALLchr.txt.gz", head=FALSE, stringsAsFactors=FALSE))
DATA_72 = as.data.frame(data.table::fread("/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/EGAD00001002672_ALL.ALLchr.txt.gz", head=FALSE, stringsAsFactors=FALSE))
DATA_73 = as.data.frame(data.table::fread("/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS/EGAD00001002673_ALL.ALLchr.txt.gz", head=FALSE, stringsAsFactors=FALSE))

DATA=DATA_73
colnames(DATA) = c("idx1","chr1","start1","end1","id1","idx2","chr2","start2","end2","id2","corr","pval")
Q = qvalue(DATA$pval)
DATA$qval = Q$qvalues
DATAs = DATA[Q$qvalue < 0.01, ]


LINKS.TRH = data.frame(from=DATAs$id1, to=DATAs$id2, weigth=ifelse(DATAs$corr>0, 1, 2), stringsAsFactors=FALSE)
NODES.TRH = data.frame(id=names(table(c(LINKS.TRH$from, LINKS.TRH$to))), chr=matrix(unlist(strsplit(names(table(c(LINKS.TRH$from, LINKS.TRH$to))), split="_")), ncol=3, byrow=TRUE)[, 1], stringsAsFactors=FALSE)
DATAnet.TRH = graph_from_data_frame(d=LINKS.TRH, vertices=NODES.TRH, directed=FALSE)
communities  = fastgreedy.community(DATAnet.TRH)
N_TRH = max(communities$membership)
N_TRH_TOSHOW = 50

###### PLOT: CALL TRHs using greedy algorithm
pdf("TRH_group_sizes_72.pdf")

df = data.frame(ID=1:N_TRH,SIZE=rle(sort(communities$membership))$length,stringsAsFactors=F)
df <-df[order(-df$SIZE),]
df <- data.frame(head(df$SIZE,50),seq(1,50))
colnames(df)=c("SIZE","ID")
p2<- ggplot(df, aes(x = ID, y = SIZE)) + geom_bar(stat = "identity", fill="steelblue", width=0.7) + labs(title=paste0(N_TRH_TOSHOW," TRHs out of ",N_TRH), x ="TRH index", y = "TRH Size") 
p2 + theme(axis.text=element_text(size=16), axis.title=element_text(size=16) ) + theme_minimal() + scale_y_continuous(trans = 'log10')


pdf("TRH_group_sizes_73.pdf") # max is 31

df = data.frame(ID=1:N_TRH,SIZE=rle(sort(communities$membership))$length,stringsAsFactors=F)
for (i in length(df$ID)+1:50) {
  print(i)
  df <- df %>% add_row(ID = i, SIZE = 0)
}


df <-df[order(-df$SIZE),]
df <- data.frame(head(df$SIZE,50),seq(1,50))
colnames(df)=c("SIZE","ID")
p2<- ggplot(df, aes(x = ID, y = SIZE)) + geom_bar(stat = "identity", fill="steelblue", width=0.7) + labs(title=paste0(N_TRH_TOSHOW," TRHs out of ",N_TRH), x ="TRH index", y = "TRH Size") 
p2 + theme(axis.text=element_text(size=16), axis.title=element_text(size=16) ) + theme_minimal() + scale_y_continuous(trans = 'log10')


print(p2)
dev.off()

head(df,10)


