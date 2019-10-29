#GENOTYPE RELATED COVARIATES
SEX=read.table("/data/unige/funpopgen/odelanea/SGX/V2/data_genotypes/cov/gencord_kgpgen_kgpseq.gender", header=FALSE, stringsAsFactors=FALSE)
GEN=read.table("/data/unige/funpopgen/odelanea/SGX/V2/data_genotypes/cov/gencord_kgpgen_kgpseq.genotyped", header=FALSE, stringsAsFactors=FALSE)
SIN=read.table("/data/unige/funpopgen/odelanea/SGX/V2/data_genotypes/cov/gencord_kgpgen_kgpseq.sinergia", header=FALSE, stringsAsFactors=FALSE)
PCA=read.table("/data/unige/funpopgen/odelanea/SGX/V2/data_genotypes/cov/gencord_kgpgen_kgpseq.pca", header=TRUE, stringsAsFactors=FALSE)
colnames(SEX)=c("IND", "SEX")
colnames(GEN)=c("IND", "GEN")
colnames(SIN)=c("IND", "SIN")
PCA=data.frame(t(PCA[1:3, -1]))
colnames(PCA)=c("GPC1", "GPC2", "GPC3")
PCA$IND=rownames(PCA)
GCOV=merge(SEX, GEN, by="IND")
GCOV=merge(GCOV, SIN, by="IND")
GCOV=merge(GCOV, PCA, by="IND")

outfile=paste("/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/BLUEPRINT_SGX_94SAMPLES/SAMPLE_ID/LCL.cov", sep="")
write.table(t(GCOV), file=outfile, quote=FALSE, row.names=TRUE, col.names=FALSE,sep='\t')

