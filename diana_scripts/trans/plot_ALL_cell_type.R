library(qvalue)
library(ggplot2)
library(gplots)
library(data.table)
library(RColorBrewer)
library(reshape2)

getchromFreq <- function(){
    DATAsig = DATA[DATA$pval<1e-06,]
    chromFreq = matrix(0,22,22)
    for(i in 1:21){
        for(j in (i+1):22){
            frac_signif = sum(DATAsig$chr1 == i & DATAsig$chr2 == j)/sum(DATA$chr1 == i & DATA$chr2 == j)
            chromFreq[i,j] = frac_signif
            chromFreq[j,i] = frac_signif
        }
    }
    chromFreq
}

load('TRANS_analysis.rda')
myhist_bg.NEU = hist(PCHiC$Neu,breaks = c(0,10,20,50,100,2000),plot=F)
myhist_signif.NEU = hist(PCHiC$Neu[hic_validated<0.05],breaks = c(0,10,20,50,100,2000),plot=F)
DATAs.NEU = DATA[Q$qvalue < 0.01,]
chromFreq.NEU = getchromFreq()

load('../../EGAD00001002672_CLOMICS_v3.0/TRANS/TRANS_analysis.rda')
myhist_bg.MON = hist(PCHiC$Mon,breaks = c(0,10,20,50,100,2000),plot=F)
myhist_signif.MON = hist(PCHiC$Mon[hic_validated<0.05],breaks = c(0,10,20,50,100,2000),plot=F)
DATAs.MON = DATA[Q$qvalue < 0.01,]
chromFreq.MON = getchromFreq()


load('../../EGAD00001002673_CLOMICS_v3.0/TRANS/TRANS_analysis.rda')
myhist_bg.TCL = hist(PCHiC$nCD4,breaks = c(0,10,20,50,100,2000),plot=F)
myhist_signif.TCL = hist(PCHiC$nCD4[hic_validated<0.05],breaks = c(0,10,20,50,100,2000),plot=F)
DATAs.TCL = DATA[Q$qvalue < 0.01,]
chromFreq.TCL = getchromFreq()

stop()

pdf("HiC_validation_ALL_CellTypes.pdf",5,7)
toplot = data.frame(Number = c("0-10","11-20","21-50","51-100",">100"),NEU = myhist_signif.NEU$counts/myhist_bg.NEU$counts*100, MON = myhist_signif.MON$counts/myhist_bg.MON$counts*100,
TCL = myhist_signif.TCL$counts/myhist_bg.TCL$counts*100)
toplot.melted = reshape2::melt(toplot,id.vars = "Number", measure.vars = c("NEU","MON","TCL"))
colnames(toplot.melted)[2] = "CellType"
toplot.melted$Number = factor(toplot$Number,levels = c("0-10","11-20","21-50","51-100",">100"))
g <- ggplot(toplot.melted, aes(x = Number, y = value,fill=CellType))+ ggtitle("") +
  geom_bar(position="dodge", stat="identity") + labs(x = "PC HiC Score",y = "Fraction of associations (%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))
print(g)
dev.off()

pdf("Connectivity_ALL_CellTypes.pdf", 6, 6)
hist.NEU = hist(table(c(DATAs.NEU$id1, DATAs.NEU$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
hist.MON = hist(table(c(DATAs.MON$id1, DATAs.MON$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)
hist.TCL = hist(table(c(DATAs.TCL$id1, DATAs.TCL$id2)), breaks=c(0,5,10,20,50,100,200,500),plot=F)

frac.NEU = hist.NEU$counts/sum(hist.NEU$counts)*100
frac.MON = hist.MON$counts/sum(hist.MON$counts)*100
frac.TCL = hist.TCL$counts/sum(hist.TCL$counts)*100

toplot = data.frame(NEU = frac.NEU,MON = frac.MON,TCL = frac.TCL,Number = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
toplot.melted = reshape2::melt(toplot,id.vars = "Number", measure.vars = c("NEU","MON","TCL"))
colnames(toplot.melted)[2] = "CellType"
toplot.melted$Number = factor(toplot$Number,levels = c("1-5","6-10","11-20","21-50","51-100","101-200","200-500"))
g <- ggplot(toplot.melted, aes(x = Number, y = value,fill=CellType))+ ggtitle("") +
  geom_bar(position="dodge", stat="identity") +
  labs(x = "Connectivity",y = "Fraction of CRDs (%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text.x = element_text(size = 20, angle = 45, hjust = 1),axis.text.y = element_text(size = 20))
 print(g)
dev.off()


chromFreq = matrix(0,68,68)
chromFreq.NEU.raw = chromFreq.NEU
chromFreq.MON.raw = chromFreq.MON
chromFreq.TCL.raw = chromFreq.TCL
chromFreq.NEU = chromFreq.NEU/mean(chromFreq.NEU)
chromFreq.MON = chromFreq.MON/mean(chromFreq.MON)
chromFreq.TCL = chromFreq.TCL/mean(chromFreq.TCL)
chromFreq[1:22,1:22] = log2(chromFreq.NEU)
chromFreq[24:45,24:45] = log2(chromFreq.MON)
chromFreq[47:68,47:68] = log2(chromFreq.TCL)
chromFreq[24:45,1:22] = log2(chromFreq.MON/chromFreq.NEU)
chromFreq[1:22,24:45] = log2(chromFreq.NEU/chromFreq.MON)
chromFreq[47:68,1:22] = log2(chromFreq.TCL/chromFreq.NEU)
chromFreq[1:22,47:68] = log2(chromFreq.NEU/chromFreq.TCL)
chromFreq[47:68,24:45] = log2(chromFreq.TCL/chromFreq.MON)
chromFreq[24:45,47:68] = log2(chromFreq.MON/chromFreq.TCL)

chromFreq[is.infinite(chromFreq)] = 0

coul <- colorRampPalette(brewer.pal(9, "RdBu"))(15)


pdf("Frequency_by_chrom_pairs_ALL.pdf")
heatmap.2(chromFreq,dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',col=coul,labRow = FALSE, labCol = FALSE,density.info='none')
dev.off()
