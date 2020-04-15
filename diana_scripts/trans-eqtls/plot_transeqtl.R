library(qvalue)

importfile <- function(CellType,eQTLmode){
    tmp = read.table(paste0(eQTLmode,'/transeqtl_',CellType,'.out'))
    Q = qvalue(tmp$V10);
    tmp$qval = Q$qvalue
    tmp
}

CellTypes = c("NEU","MON","TCEL")
CellTypes.names = c("NEU","MON","TCL")
eQTLmodes = c("aCRD","sCRD","ciseQTL")

for(mymode in eQTLmodes){
    NEU = importfile(CellTypes[1],mymode)
    MON = importfile(CellTypes[2],mymode)
    TCL = importfile(CellTypes[3],mymode)
    pdf(paste0("Boxplot_",mymode,".pdf"),height=7,width=5)
    par(mar=c(5.1, 5.1, 4.1, 2.1))
    boxplot(-log10(NEU$V10),-log10(MON$V10),-log10(TCL$V10),col=c("orchid","lightblue","palegreen3"),notch=T,names=CellTypes,ylab=c("-log10(p-value)"),cex.axis=2,cex.lab=2, frame=F)
    points(rep(1,sum(NEU$qval<0.05)),-log10(NEU$V10[which(NEU$qval<0.05)]),col="royalblue")
    points(rep(2,sum(MON$qval<0.05)),-log10(MON$V10[which(MON$qval<0.05)]),col="royalblue")
    points(rep(3,sum(TCL$qval<0.05)),-log10(TCL$V10[which(TCL$qval<0.05)]),col="royalblue")
    dev.off()
}

NEU.aCRD = importfile(CellTypes[1],eQTLmodes[1])
MON.aCRD = importfile(CellTypes[2],eQTLmodes[1])
TCL.aCRD = importfile(CellTypes[3],eQTLmodes[1])


NEU.sCRD = importfile(CellTypes[1],eQTLmodes[2])
MON.sCRD = importfile(CellTypes[2],eQTLmodes[2])
TCL.sCRD = importfile(CellTypes[3],eQTLmodes[2])


NEU.cis = importfile(CellTypes[1],eQTLmodes[3])
MON.cis = importfile(CellTypes[2],eQTLmodes[3])
TCL.cis = importfile(CellTypes[3],eQTLmodes[3])

pdf(paste0("Boxplot_ALL.pdf"),height=5,width=13)
par(mar=c(5.1, 5.1, 4.1, 2.1))
toplot = list(-log10(NEU.aCRD$V10),-log10(MON.aCRD$V10),-log10(TCL.aCRD$V10),-log10(NEU.sCRD$V10),-log10(MON.sCRD$V10),-log10(TCL.sCRD$V10),-log10(NEU.cis$V10),-log10(MON.cis$V10),-log10(TCL.cis$V10))
boxplot(toplot,col=rep(c("orchid","lightblue","palegreen3"),3),notch=T,names=rep(CellTypes.names,3),ylab=c("-log10(p-value)"),cex.axis=2,cex.lab=2, frame=F)

points(rep(1,sum(NEU.aCRD$qval<0.05)),-log10(NEU.aCRD$V10[which(NEU.aCRD$qval<0.05)]),col="royalblue")
points(rep(2,sum(MON.aCRD$qval<0.05)),-log10(MON.aCRD$V10[which(MON.aCRD$qval<0.05)]),col="royalblue")
points(rep(3,sum(TCL.aCRD$qval<0.05)),-log10(TCL.aCRD$V10[which(TCL.aCRD$qval<0.05)]),col="royalblue")
abline(v=3.5,lty=2)

points(rep(4,sum(NEU.sCRD$qval<0.05)),-log10(NEU.sCRD$V10[which(NEU.sCRD$qval<0.05)]),col="royalblue")
points(rep(5,sum(MON.sCRD$qval<0.05)),-log10(MON.sCRD$V10[which(MON.sCRD$qval<0.05)]),col="royalblue")
points(rep(6,sum(TCL.sCRD$qval<0.05)),-log10(TCL.sCRD$V10[which(TCL.sCRD$qval<0.05)]),col="royalblue")
abline(v=6.5,lty=2)

points(rep(7,sum(NEU.cis$qval<0.05)),-log10(NEU.cis$V10[which(NEU.cis$qval<0.05)]),col="royalblue")
points(rep(8,sum(MON.cis$qval<0.05)),-log10(MON.cis$V10[which(MON.cis$qval<0.05)]),col="royalblue")
points(rep(9,sum(TCL.cis$qval<0.05)),-log10(TCL.cis$V10[which(TCL.cis$qval<0.05)]),col="royalblue")


dev.off()
