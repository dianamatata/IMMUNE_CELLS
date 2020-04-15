library(qvalue)
library(RColorBrewer)

importfile <- function(path){
    tmp = read.table(path)
    Q = qvalue(tmp$V10);
    tmp$qval = Q$qvalue
    tmp
}

CellTypes = c("NEU","MON","TCEL")
CellTypes.names = c("NEU","MON","TCL")
eQTLmodes = c("aCRD","sCRD","eGene")


NEU.aCRD = importfile("./aCRD/transeqtl_NEU.out")
NEU.sCRD = importfile("./sCRD/transeqtl_NEU.out")
NEU.cis = importfile("./ciseQTL/transeqtl_NEU.out")

MON.aCRD = importfile("../../../EGAD00001002672_CLOMICS_v3.0/TRANS-EQTL/Q0.05/aCRD/transeqtl_MON.out")
MON.sCRD = importfile("../../../EGAD00001002672_CLOMICS_v3.0/TRANS-EQTL/Q0.05/sCRD/transeqtl_MON.out")
MON.cis = importfile("../../../EGAD00001002672_CLOMICS_v3.0/TRANS-EQTL/Q0.05/ciseQTL/transeqtl_MON.out")

TCL.aCRD = importfile("../../../EGAD00001002673_CLOMICS_v3.0/TRANS-EQTL/Q0.05/aCRD/transeqtl_TCEL.out")
TCL.sCRD = importfile("../../../EGAD00001002673_CLOMICS_v3.0/TRANS-EQTL/Q0.05/sCRD/transeqtl_TCEL.out")
TCL.cis = importfile("../../../EGAD00001002673_CLOMICS_v3.0/TRANS-EQTL/Q0.05/ciseQTL/transeqtl_TCEL.out")


mycolors = colorRampPalette(brewer.pal(9,"Set1"))(9)
pdf("QQplot_ALL.pdf")
plot(sort(-log10(runif(nrow(NEU.aCRD)))), sort(-log10(NEU.aCRD$V10)), xlab="Expected -log10(P-values)", ylab="Observed -log10(P-values)",xlim=c(0,6),ylim=c(0,25),col=mycolors[1],type='b',pch = 19)
lines(sort(-log10(runif(nrow(MON.aCRD)))), sort(-log10(MON.aCRD$V10)),col=mycolors[2],type='b',pch = 19)
lines(sort(-log10(runif(nrow(TCL.aCRD)))), sort(-log10(TCL.aCRD$V10)),col=mycolors[3],type='b',pch = 19)

lines(sort(-log10(runif(nrow(MON.aCRD)))), sort(-log10(MON.aCRD$V10)),col=mycolors[4],type='b',pch = 19)
lines(sort(-log10(runif(nrow(MON.sCRD)))), sort(-log10(MON.sCRD$V10)),col=mycolors[5],type='b',pch = 19)
lines(sort(-log10(runif(nrow(MON.cis)))), sort(-log10(MON.cis$V10)),col=mycolors[6],type='b',pch = 19)

lines(sort(-log10(runif(nrow(TCL.aCRD)))), sort(-log10(TCL.aCRD$V10)),col=mycolors[7],type='b',pch = 19)
lines(sort(-log10(runif(nrow(TCL.sCRD)))), sort(-log10(TCL.sCRD$V10)),col=mycolors[8],type='b',pch = 19)
lines(sort(-log10(runif(nrow(TCL.cis)))), sort(-log10(TCL.cis$V10)),col=mycolors[9],type='b',pch = 19)

legend(x='topleft',legend=c("NEU aCRD","NEU sCRD","NEU eGene","MON aCRD","MON sCRD","MON eGene","TCL aCRD","TCL sCRD","TCL eGene"),fill = mycolors)
abline(0, 1, col="grey")
dev.off()

pdf("QQplot_aCRD.pdf")
par(mar=c(6.1, 5.1, 4.1, 2.1))
plot(sort(-log10(runif(nrow(NEU.aCRD)))), sort(-log10(NEU.aCRD$V10)), xlab="Expected -log10(P-values)", ylab="Observed -log10(P-values)",xlim=c(0,6),ylim=c(0,25),col=mycolors[1],type='o',pch = 19,cex.axis=2,cex.lab=2)
lines(sort(-log10(runif(nrow(MON.aCRD)))), sort(-log10(MON.aCRD$V10)),col=mycolors[2],type='o',pch = 19)
lines(sort(-log10(runif(nrow(TCL.aCRD)))), sort(-log10(TCL.aCRD$V10)),col=mycolors[3],type='o',pch = 19)
legend(x='topleft',legend=c("NEU","MON","TCL"),fill = mycolors[1:3],cex=1.5)
abline(0, 1, col="grey")
dev.off()

pdf("QQplot_sCRD.pdf")
par(mar=c(6.1, 5.1, 4.1, 2.1))
plot(sort(-log10(runif(nrow(NEU.sCRD)))), sort(-log10(NEU.sCRD$V10)), xlab="Expected -log10(P-values)", ylab="Observed -log10(P-values)",xlim=c(0,4.8),col=mycolors[1],type='o',pch = 19,cex.axis=2,cex.lab=2)
lines(sort(-log10(runif(nrow(MON.sCRD)))), sort(-log10(MON.sCRD$V10)),col=mycolors[2],type='o',pch = 19)
lines(sort(-log10(runif(nrow(TCL.sCRD)))), sort(-log10(TCL.sCRD$V10)),col=mycolors[3],type='o',pch = 19)
legend(x='topleft',legend=c("NEU","MON","TCL"),fill = mycolors[1:3],cex=1.5)
abline(0, 1, col="grey")
dev.off()

pdf("QQplot_eGene.pdf")
par(mar=c(6.1, 5.1, 4.1, 2.1))
plot(sort(-log10(runif(nrow(NEU.cis)))), sort(-log10(NEU.cis$V10)), xlab="Expected -log10(P-values)", ylab="Observed -log10(P-values)",xlim=c(0,7),ylim=c(0,25),col=mycolors[1],type='o',pch = 19,cex.axis=2,cex.lab=2)
lines(sort(-log10(runif(nrow(MON.cis)))), sort(-log10(MON.cis$V10)),col=mycolors[2],type='o',pch = 19)
lines(sort(-log10(runif(nrow(TCL.cis)))), sort(-log10(TCL.cis$V10)),col=mycolors[3],type='o',pch = 19)
legend(x='topleft',legend=c("NEU","MON","TCL"),fill = mycolors[1:3],cex=1.5)
abline(0, 1, col="grey")
dev.off()


pdf(paste0("Boxplot_ALL_ALL.pdf"),height=8,width=15)
par(mar=c(5.1, 5.1, 4.1, 2.1))
toplot = list(-log10(NEU.aCRD$V10),-log10(NEU.sCRD$V10),-log10(NEU.cis$V10),-log10(MON.aCRD$V10),-log10(MON.sCRD$V10),-log10(MON.cis$V10),-log10(TCL.aCRD$V10),-log10(TCL.sCRD$V10),-log10(TCL.cis$V10))
boxplot(toplot,col=rep(c("orchid","lightblue","palegreen3"),each=3,times=1),notch=T,names=rep(eQTLmodes,times=3),ylab=c("-log10(p-value)"),cex.axis=2,cex.lab=2, frame=F)

points(rep(1,sum(NEU.aCRD$qval<0.05)),-log10(NEU.aCRD$V10[which(NEU.aCRD$qval<0.05)]),col="royalblue")
points(rep(2,sum(NEU.sCRD$qval<0.05)),-log10(NEU.sCRD$V10[which(NEU.sCRD$qval<0.05)]),col="royalblue")
points(rep(3,sum(NEU.cis$qval<0.05)),-log10(NEU.cis$V10[which(NEU.cis$qval<0.05)]),col="royalblue")
abline(v=3.5,lty=2)

points(rep(4,sum(MON.aCRD$qval<0.05)),-log10(MON.aCRD$V10[which(MON.aCRD$qval<0.05)]),col="royalblue")
points(rep(5,sum(MON.sCRD$qval<0.05)),-log10(MON.sCRD$V10[which(MON.sCRD$qval<0.05)]),col="royalblue")
points(rep(6,sum(MON.cis$qval<0.05)),-log10(MON.cis$V10[which(MON.cis$qval<0.05)]),col="royalblue")
abline(v=6.5,lty=2)

points(rep(7,sum(TCL.aCRD$qval<0.05)),-log10(TCL.aCRD$V10[which(TCL.aCRD$qval<0.05)]),col="royalblue")
points(rep(8,sum(TCL.sCRD$qval<0.05)),-log10(TCL.sCRD$V10[which(TCL.sCRD$qval<0.05)]),col="royalblue")
points(rep(9,sum(TCL.cis$qval<0.05)),-log10(TCL.cis$V10[which(TCL.cis$qval<0.05)]),col="royalblue")


dev.off()
