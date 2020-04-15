library(ggplot2)
library(reshape2)

load('Gene_CRD_associations.RDAta')
mat_hic.NEU = mat_hic
coexpressed_mat_hic.NEU = coexpressed_mat_hic

load('../../EGAD00001002672_CLOMICS_v3.0/mapping_aCRD_gene/Gene_CRD_associations.RDAta') # does not exist anymore,  find . -name "*RData"
mat_hic.MON = mat_hic
coexpressed_mat_hic.MON = coexpressed_mat_hic

load('../../EGAD00001002673_CLOMICS_v3.0/mapping_aCRD_gene/Gene_CRD_associations.RDAta')
mat_hic.TCL = mat_hic
coexpressed_mat_hic.TCL = coexpressed_mat_hic

pdf("HiC_support_gene_CRD_associations_ALL_CellTypes.pdf")
toplot = data.frame(NEU = mat_hic.NEU[,1],MON=mat_hic.MON[,1],TCL=mat_hic.TCL[,1],dist = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
toplot.melted = reshape2::melt(toplot,id.vars = "dist", measure.vars = c("NEU","MON","TCL"))
colnames(toplot.melted)[2] = "CellType"
toplot.melted$dist <- factor(toplot$dist ,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))

ggplot(toplot.melted, aes(x = dist, y = value,group=CellType))+ ggtitle("HiC support for gene-CRD associations") +
geom_line(aes(color=CellType),size=1) + geom_point(aes(color=CellType)) +
labs(x = "Distance",y = "Fraction with HiC support (%)") +
   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

pdf("HiC_support_gene_pairs_with_same_CRD_ALL_CellTypes.pdf")
toplot = data.frame(NEU = coexpressed_mat_hic.NEU[,1],MON = coexpressed_mat_hic.MON[,1],TCL = coexpressed_mat_hic.TCL[,1],dist = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))
toplot.melted = reshape2::melt(toplot,id.vars = "dist", measure.vars = c("NEU","MON","TCL"))
colnames(toplot.melted)[2] = "CellType"
toplot.melted$dist <- factor(toplot$dist ,levels = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))

ggplot(toplot.melted, aes(x = dist, y = value,fill=CellType))+ ggtitle("") +
   geom_bar(position="dodge", stat="identity") +
  labs(x = "Distance between co-expressed genes",y = "Fraction with HiC support (%)") +
   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
