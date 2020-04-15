# mozaic plot for Fig 6.5: Overlap between trans eGenes (aCRD, sCRD and eGene) across the three immune cell types.

library(UpSetR)
library(reshape2)
library(ggplot2)
library(RColorBrewer)


NEU = scan('/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS-EQTL/Q0.05/70_trans_eGenes.txt',what='')
MON = scan('/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS-EQTL/Q0.05/72_trans_eGenes.txt',what='')
TCL = scan('/Users/dianaavalos/Programming/THREE_CELL_TYPES__CLOMICS__EGAD00001002670_CLOMICS_v3.0__TRANS-EQTL/Q0.05/73_trans_eGenes.txt',what='')

N_M = length(intersect(NEU,MON))
N_T = length(intersect(NEU,TCL))
M_T = length(intersect(TCL,MON))
N_M_T = length(intersect(NEU,intersect(TCL,MON)))

expressionInput <- c(NEU = length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T), MON = length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T), TCL = length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T), `NEU&MON` = N_M-N_M_T,
 `NEU&TCL` = N_T-N_M_T,`MON&TCL` = M_T-N_M_T, `NEU&MON&TCL` = N_M_T)

pdf("Overlap_gene_sets.pdf",paper='a4r')
upset(fromExpression(expressionInput), order.by = "degree",text.scale=1.8)
upset(fromExpression(expressionInput), order.by = "freq",text.scale=1.8,main.bar.color="orange",matrix.color="blue")
dev.off()

toplot = data.frame(NEU=c(length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T),0,0,N_M-N_M_T,N_T-N_M_T,0,N_M_T),MON=c(0,length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T),0,N_M-N_M_T,0,M_T-N_M_T,N_M_T),
TCL=c(0,0,length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T),0,N_T-N_M_T,M_T-N_M_T,N_M_T),Group =c("NEU","MON","TCL","NEU+MON","NEU+TCL","MON+TCL","ALL"))

toplot.molten = melt(toplot,id.vars="Group")
toplot.molten$Group <- factor(toplot.molten$Group ,levels = c("NEU","MON","TCL","NEU+MON","NEU+TCL","MON+TCL","ALL"))

pdf("Overlap_stacked_barplot.pdf")
g <-ggplot(toplot.molten, aes(fill=Group, y=value, x=variable)) + geom_bar(positio="stack", stat="identity") + labs(x = "Cell Type",y = "Trans eGenes") + scale_fill_brewer(palette = "Set3") +
theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))
print(g)
dev.off()


# Diana 
NEUT = length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T)
MONO = length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T)
TCEL = length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T)
`NEUT&MONO` = N_M-N_M_T
`NEUT&TCEL` = N_T-N_M_T
`MONO&TCEL` = M_T-N_M_T
`NEUT&MONO&TCEL` = N_M_T

One_cell=NEUT+MONO+TCEL
Two_cells=`NEUT&MONO`+`NEUT&TCEL`+`MONO&TCEL`
Three_cells=NEUT&MONO&TCEL

expressionInput2 <- c(One_cell=NEUT+MONO+TCEL, Two_cells=`NEUT&MONO`+`NEUT&TCEL`+`MONO&TCEL`,Three_cells=NEUT&MONO&TCEL )

toplot2 = data.frame(NEU=c(length(NEU)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T),N_M-N_M_T+N_T-N_M_T+0,N_M_T),MON=c(length(MON)-N_M_T-(N_M-N_M_T)-(N_T-N_M_T),N_M-N_M_T+0+M_T-N_M_T,N_M_T),
                    TCL=c(length(TCL)-N_M_T-(N_T-N_M_T)-(M_T-N_M_T),0+N_T-N_M_T+M_T-N_M_T,N_M_T),Group =c("1cell_type","2cell_types","3cell_types"))

## plots
toplot2.molten = melt(toplot2,id.vars="Group")
toplot2.molten$Group <- factor(toplot2.molten$Group ,levels = c("1cell_type","2cell_types","3cell_types"))
g <-ggplot(toplot2.molten, aes(fill=Group, y=value, x=variable)) + geom_bar(positio="stack", stat="identity") + labs(x = "Cell Type",y = "Trans eGenes") + scale_fill_brewer(palette = "Set3") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))
print(g)

## plots
#optional to normalize
topoplot_normalized=apply(toplot2[1:3],2,norm<-function(x){return (x/sum(x)*100)})
topoplot_normalized <- as.data.frame((topoplot_normalized))
topoplot_normalized$Group<-c("1cell_type","2cell_types","3cell_types")

toplot3.molten = melt(topoplot_normalized,id.vars="Group")
toplot3.molten$Group <- factor(toplot3.molten$Group ,levels = c("1cell_type","2cell_types","3cell_types"))
g <-ggplot(toplot3.molten, aes(fill=Group, y=value, x=variable)) + geom_bar(positio="stack", stat="identity") + labs(x = "Cell Type",y = "Trans eGenes") + scale_fill_brewer(palette = "Set3") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))
print(g)


## mozaic plot not ready yet
library(ggmosaic)
mosaic(HairEyeColor, shade=TRUE, legend=TRUE)
g <-ggplot(data = toplot2.molten)
print(g)
+ geom_mosaic(aes(x = product(variable, value), fill=value), na.rm=TRUE) + 
  labs(x = "Is it rude recline? ", title='f(DoYouRecline | RudeToRecline) f(RudeToRecline)')

ggplot(data = fly) +
  geom_mosaic(aes(x = product(DoYouRecline, RudeToRecline), fill=DoYouRecline), na.rm=TRUE) + 
  labs(x = "Is it rude recline? ", title='f(DoYouRecline | RudeToRecline) f(RudeToRecline)')
