library(ggplot2)
library(reshape2)

# path:/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0/mapping_aCRD_gene%


load('/Users/dianaavalos/Programming/Hi-C_correlated_peaks_70/Correlation_vs_PCHiC_70.RData')
myhist_bg.NEU = hist(PCHiC$Neu[hic_validated<1],breaks = c(seq(0,20,2),2000),plot=F)
myhist_signif.NEU = hist(PCHiC$Neu[hic_validated<0.01],breaks = c(seq(0,20,2),2000),plot=F)

load('/Users/dianaavalos/Programming/Hi-C_correlated_peaks_70/Correlation_vs_PCHiC_72.RData')
myhist_bg.MON = hist(PCHiC$Mon[hic_validated<1],breaks = c(seq(0,20,2),2000),plot=F)
myhist_signif.MON = hist(PCHiC$Mon[hic_validated<0.01],breaks = c(seq(0,20,2),2000),plot=F)

load('/Users/dianaavalos/Programming/Hi-C_correlated_peaks_70/Correlation_vs_PCHiC_73.RData')
myhist_bg.TCL = hist(PCHiC$nCD4[hic_validated<1],breaks = c(seq(0,20,2),2000),plot=F)
myhist_signif.TCL = hist(PCHiC$nCD4[hic_validated<0.01],breaks = c(seq(0,20,2),2000),plot=F)

toplot = data.frame(NEU = myhist_signif.NEU$counts/myhist_bg.NEU$counts*100,MON=myhist_signif.MON$counts/myhist_bg.MON$counts*100,
    TCL=myhist_signif.TCL$counts/myhist_bg.TCL$counts*100,Number = c(seq(2,20,2),">20"))

toplot.molten = melt(toplot,id.vars="Number")
toplot.molten$Number = factor(toplot.molten$Number,levels = c(seq(2,20,2),">20"))
colnames(toplot.molten)[2] = "CellType"

pdf("HiC_validation_ALL_CellTypes.pdf",7,5)
g <- ggplot(toplot.molten, aes(x = Number, y = value,fill=CellType))+ ggtitle("") +
   geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =  c("#F1BB7B", "#FD6467", "#5B1A18")) +
   theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
   theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 0, hjust = 1))+
   labs(x = "PCHiC Score",y = "Fraction correlated peaks (%)") 
print(g)
dev.off()

pdf("HiC_validation_ALL_CellTypes_linePlot.pdf",7,5)
g <- ggplot(toplot.molten, aes(x = Number, y = value,group=CellType))+ ggtitle("") +
   geom_line(aes(color=CellType)) + geom_point(aes(color=CellType))+
   scale_colour_manual(values=c(NEU="#F1BB7B",MON="#FD6467",TCL="#5B1A18"))+
 #  geom_line(aes(color=CellType)) + geom_point(aes(color=CellType)) +
  labs(x = "PCHiC Score",y = "Fraction correlated peaks (%)") +
   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 0, hjust = 1))
   print(g)
dev.off()
