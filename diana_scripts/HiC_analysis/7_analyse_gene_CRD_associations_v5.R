library(qvalue)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(tidyverse)

# path:/srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0/mapping_aCRD_gene%
# /srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CORRELATION/expression

allCRDs = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks_70/EGAD00001002670_ALL.modules.MOD1.NRE2.txt.gz',header=F)
protein_coding_genes = scan("/Users/dianaavalos/Programming/reference_files/gencode.v15.annotation.protein_coding.gene_id.txt",what="")
long_nc_RNA_genes = scan("/Users/dianaavalos/Programming/reference_files/gencode.v15.annotation.long_noncoding_RNAs.gene_id.txt",what="")
genelist = c(protein_coding_genes,long_nc_RNA_genes)
mapdata = read.table('/Users/dianaavalos/Programming/Hi-C_correlated_peaks_70/mapping_gene_CRD_mean_ALL_70.txt',stringsAsFactors=F)

colnames(mapdata) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
 "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","nominal_pval","slope","top_variant")
#Filter protein coding and long nc RNA genes
mapdata = mapdata[mapdata$phenotype_ID %in% genelist,]
mapdata = mapdata[order(mapdata[,2],mapdata[,3]),]

nb_CRD_not_associated = nrow(allCRDs) - length(unique(mapdata$CRD_ID))

d = read.table("/Users/dianaavalos/Programming/Hi-C_correlated_peaks_70/mapping_gene_CRD_mean_ALL_70.txt", hea=F, stringsAsFactors=F)

colnames(d) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
"nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","degree_freedom",
"Dummy","1st_param_beta","2nd_param_beta","nominal_pval","slope","empirical_pval","beta_pval")

# attribut 'names' [19] doit être de même longueur que le vecteur [14]. next one of 14, not sure which column names to keep though
colnames(d) = c("phenotype_ID","phenotype_ID_chr","phenotype_ID_start","phenotype_ID_end","phenotype_ID_strand",
                "nb_variants","distance","CRD_ID","CRD_ID_chr","CRD_ID_start","CRD_ID_end","degree_freedom",
                "Dummy","1st_param_beta")

d = d[d$phenotype_ID %in% genelist,]

cordata = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks_70/EGAD00001002675_RNA.ALL.txt.gz')
cordata$distance = cordata$V4-cordata$V3
cordata = cordata[,c(3,4,5,6,8,9)]
colnames(cordata)[1:5] = c("pos1","pos2","gene1","gene2","pval")
cordata$pval.adj = p.adjust(cordata$pval,method='fdr')
cordata$sameCRD = -1
setkeyv(cordata,c("gene1","gene2"))

cordata = cordata[cordata$gene1 %in% genelist,]
cordata = cordata[cordata$gene2 %in% genelist,]

cordata$rowid = c(1:nrow(cordata))


nb_genes_not_associated = length(unique(cordata$gene1)) - length(unique(mapdata$phenotype_ID))

# /data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/EXTERNAL_DATA/PCHiC_peak_matrix_cutoff5.tsv: same for all cells
PCHiC = fread('/Users/dianaavalos/Programming/Hi-C_correlated_peaks_70/PCHiC_peak_matrix_cutoff5.tsv')
colnames(PCHiC)[1] = "baitChr"

baitbed <- GRanges(seqnames=PCHiC$baitChr,ranges=IRanges(start=PCHiC$baitStart, end=PCHiC$baitEnd))
oebed <- GRanges(seqnames=PCHiC$oeChr,ranges=IRanges(start=PCHiC$oeStart, end=PCHiC$oeEnd))

genebed <- GRanges(seqnames=mapdata$phenotype_ID_chr,ranges=IRanges(start=mapdata$phenotype_ID_start, end=mapdata$phenotype_ID_end))
CRDbed <- GRanges(seqnames=mapdata$CRD_ID_chr,ranges=IRanges(start=mapdata$CRD_ID_start, end=mapdata$CRD_ID_end))

#fwd
x = findOverlaps(baitbed,genebed)
y = findOverlaps(oebed,CRDbed)
tmp = rbind(as.data.frame(x),as.data.frame(y))
validated.fwd = tmp[which(duplicated(tmp)),]

#bwd
x = findOverlaps(baitbed,CRDbed)
y = findOverlaps(oebed,genebed)
tmp = rbind(as.data.frame(x),as.data.frame(y))
validated.bwd = tmp[which(duplicated(tmp)),]

validated = unique(rbind(validated.fwd,validated.bwd))

hic_validated = rep(1,nrow(PCHiC))

for(i in 1:nrow(validated)){
    currenthic = validated$queryHits[i]
    currentpval = mapdata$nominal_pval[validated$subjectHits[i]]
    if(currentpval<hic_validated[currenthic]){
        hic_validated[currenthic] = currentpval
    }
}

mapdata_validated = rep(0,nrow(mapdata))

for(i in 1:nrow(validated)){
    currenthic = validated$queryHits[i]
    currentmap = validated$subjectHits[i]
    currentHiCScore = PCHiC$Neu[currenthic]
    if(currentHiCScore>mapdata_validated[currentmap]){
        mapdata_validated[currentmap] = currentHiCScore
    }
}

pdf("Distance_analysis.pdf")
hist(mapdata$distance[which(mapdata$distance!=0)],n=101,main='Distance between genes and CRDs',xlab='Distance [bp]',cex.lab=1.3,cex.axis=1.3)
relative_position = (mapdata$phenotype_ID_start[which(mapdata$distance==0)]-mapdata$CRD_ID_start[which(mapdata$distance==0)])/(mapdata$CRD_ID_end[which(mapdata$distance==0)]-mapdata$CRD_ID_start[which(mapdata$distance==0)])
hist(relative_position,xlim=c(0,1),n=50, main='Relative positive of genes within CRDs',xlab='Relative position',cex.lab=1.3,cex.axis=1.3)
dev.off()

pdf("Histogram_effect_size.pdf")
hist(mapdata$slope,n=101,main='Effect size of the associations between genes and CRDs',xlab='Slope',cex.lab=1.3,cex.axis=1.3)
dev.off()

pdf("Histogram_pvalue.pdf")
hist(d$beta_pval,main=paste0("Histogram of adjusted p-values%\n",nrow(mapdata)," significant associations at 5% FDR"),xlab="P-values",cex.lab=1.3,cex.axis=1.3)
dev.off()

pdf("Boxplot_pvalue_by_distance.pdf")
within_CRD = -log10(mapdata$nominal_pval[which(mapdata$distance==0)])
less_1kb = -log10(mapdata$nominal_pval[which(mapdata$distance>0 & mapdata$distance<1e03)])
less_10kb = -log10(mapdata$nominal_pval[which(mapdata$distance>1e03 & mapdata$distance<1e04)])
less_100kb = -log10(mapdata$nominal_pval[which(mapdata$distance>1e04 & mapdata$distance<1e05)])
less_1Mb = -log10(mapdata$nominal_pval[which(mapdata$distance>1e05 & mapdata$distance<1e06)])
boxplot(within_CRD,less_1kb,less_10kb,less_100kb,less_1Mb,names=c("within_CRD","<1kb","<10kb","<100kb","<1Mb"),ylab="-log10(p-value)",cex.lab=1.3,cex.axis=1.3)
dev.off()
#stop()

# Fig 3.6:  Number of genes or CRDs as a function of the number of CRDs and vice versa, respectively. - in neutrophils

pdf("Connectivity_CRD_gene.pdf")
CRDhist = hist(table(mapdata$CRD_ID),breaks=c(0,1,2,3,4,100),plot=F)
genehist = hist(table(mapdata$phenotype_ID),breaks=c(0,1,2,3,100),plot=F)
barplot(c(nb_CRD_not_associated,CRDhist$counts),names=c("0","1","2","3","4","5+"),main="Number of genes associated with each CRD",cex.names=1.5,cex.axis=1.5)
barplot(c(nb_genes_not_associated,genehist$counts),names=c("0","1","2","3","4+"),main="Number of CRDs associated with each gene",cex.names=1.5,cex.axis=1.5)

ggplot(data.frame(counts = c(nb_CRD_not_associated,CRDhist$counts),Number = c("0","1","2","3","4","5+")), aes(x = Number, y = counts))+ ggtitle("Genes associated with each CRD") +
  geom_bar(stat = "identity",fill="#5B1A18") +
  geom_text(aes(label = sprintf("%.2f%%", counts/sum(counts) * 100)),vjust = -.5, size =6) + labs(x = "Number of associated genes",y = "CRD counts") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))

 ggplot(data.frame(counts = c(nb_genes_not_associated,genehist$counts),Number = c("0","1","2","3","4+")), aes(x = Number, y = counts))+ ggtitle("CRDs associated with each gene") +
    geom_bar(stat = "identity",fill="#F1BB7B") +
    geom_text(aes(label = sprintf("%.2f%%", counts/sum(counts) * 100)),vjust = -.5, size =6) + labs(x = "Number of associated CRDs",y = "Gene counts") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))

dev.off()


# for computing the Hi-C support
#weird because mapdata$distance is the dist between CRDs and genes
crd_dist_bins = c(0,1,1e04,2e04,5e04,1e05,2e05,5e05,1e06)

compute_ratio_hic_mapdata <-function(CRDmindist,CRDmaxdist,cutoff=5){
    up = sum(abs(mapdata$distance)>=CRDmindist & abs(mapdata$distance)<CRDmaxdist & mapdata_validated>cutoff,na.rm=T)
    down = sum(abs(mapdata$distance)>=CRDmindist & abs(mapdata$distance)<CRDmaxdist,na.rm=T)
    up/down*100
}

mat_hic = matrix(0,nrow=(length(crd_dist_bins)-1),ncol=1)

for(i in 1:(length(crd_dist_bins)-1)){
        mat_hic[i,1] = compute_ratio_hic_mapdata(crd_dist_bins[i],crd_dist_bins[i+1])
}

###########

## here for the 4.3 plot

pdf("HiC_support_gene_CRD_associations.pdf")
toplot = data.frame(counts = mat_hic[,1],dist = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
toplot$dist <- factor(toplot$dist ,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))

ggplot(toplot, aes(x = dist, y = counts))+ ggtitle("HiC support for gene-CRD associations") +
   geom_bar(stat = "identity",fill="#56B4E9") +
   labs(x = "Distance",y = "Fraction with HiC support (%)") +
   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

#### 4.3 fig barplot
toplot2 = data.frame(NEU=c(33.01731, 26.62890, 24.71910, 21.97393, 23.14815, 22.64875, 15.50095, 11.81319),
                    MON=c(38.554217, 33.132530, 24.719101, 29.946524, 37.068966, 39.062500, 25.619835,  8.108108),
                    TCL=c(33.526235, 23.300971, 23.102310, 26.016260, 23.843416, 25.641026, 19.238901,  9.234828),
                    dist=c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))


toplot2.molten = melt(toplot2,id.vars="dist")
toplot2.molten$dist = factor(toplot2.molten$dist,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
colnames(toplot2.molten)[2] = "CellType"

g <- ggplot(toplot2.molten, aes(x = dist, y = value,fill=CellType))+ ggtitle("HiC support for gene-CRD associations") +
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =  c("#F1BB7B", "#FD6467", "#5B1A18")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Distance",y = "Fraction with HiC support (%)") 
print(g)




#######################
#Co-expression analysis
#######################

#Annotate cordata with CRD mapping

CRD_IDs = unique(mapdata$CRD_ID)
IDXs = c()
CRDdistance = c()
hicsupport = c()
for(i in 1:length(CRD_IDs)){
    # cat(i,"\n")
    idx = which(mapdata$CRD_ID == CRD_IDs[i])
    if(length(idx)>1){
        for(k in 1:(length(idx)-1)){
        	for(l in (k+1):length(idx)){
                	geneA = mapdata$phenotype_ID[idx[k]]
                	geneB = mapdata$phenotype_ID[idx[l]]
			myhit = c(unlist(cordata[gene1 == geneB & gene2 == geneA,9]),unlist(cordata[gene1 == geneA & gene2 == geneB,9]))
			if(length(myhit)>0){
                #cat(CRD_IDs[i],geneA,geneB,sep='\n')
				IDXs = c(IDXs,myhit)
				tmp = mean(abs(mapdata$distance[idx[k]]),abs(mapdata$distance[idx[l]]))
                tmp2 = mean(mapdata_validated[idx[k]],mapdata_validated[idx[l]])
				CRDdistance = c(CRDdistance,tmp)
                hicsupport = c(hicsupport,tmp2)

            }
		}
        }
    }
}

tmpdf = data.table(idx=IDXs,dist=CRDdistance,hicsupport=hicsupport)
tmpdf <- tmpdf %>%
    group_by(idx) %>%
    summarize(mindist=min(dist),hic=mean(hicsupport))

cordata$sameCRD[tmpdf$idx] = tmpdf$mindist
cordata$hic = 0
cordata$hic[tmpdf$idx] = tmpdf$hic


compute_ratio <-function(coexpressed=T,mindist,maxdist,CRDmindist,CRDmaxdist,pval.cutoff=0.01){
    if(coexpressed){
        up = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
        down = sum(cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    } else {
        up = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
        down = sum(cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    }
    up/down*100
}

compute_OR <-function(mindist,maxdist,pval.cutoff=0.01){
    coexpressed_same_CRD = sum(cordata$sameCRD>=0 & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    coexpressed_notsame_CRD = sum(!(cordata$sameCRD>=0) & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    notcoexpressed_same_CRD = sum(cordata$sameCRD>=0 & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    notcoexpressed_notsame_CRD = sum(!(cordata$sameCRD>=0) & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    m = matrix(c(coexpressed_same_CRD,coexpressed_notsame_CRD,notcoexpressed_same_CRD,notcoexpressed_notsame_CRD),ncol=2)
    pval = fisher.test(m)$p.value
    OR = fisher.test(m)$estimate
    list(OR=OR,pval=pval)
}


gene_dist_bins = c(0,1e04,2e04,5e04,1e05,2e05,5e05,1e06)
crd_dist_bins = c(0,1,1e04,1e05,1e06)
coexpressed_mat = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=(length(crd_dist_bins)-1))
notcoexpressed_mat = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=(length(crd_dist_bins)-1))

for(i in 1:(length(gene_dist_bins)-1)){
    for(k in 1:(length(crd_dist_bins)-1)){
        coexpressed_mat[i,k] = compute_ratio(T,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1])
        notcoexpressed_mat[i,k] = compute_ratio(F,gene_dist_bins[i],gene_dist_bins[i+1],crd_dist_bins[k],crd_dist_bins[k+1])
    }
}

coexpressed_mat = data.frame(coexpressed_mat)
colnames(coexpressed_mat) = c("inside","<10kb","<100kb","<1Mb")
coexpressed_mat$row = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000")
coexpressed_mat$row <- factor(coexpressed_mat$row,levels = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000"))
coexpressed_table <- melt(coexpressed_mat, id.vars = "row")
colnames(coexpressed_table) = c("Gene","CRD","value")

notcoexpressed_mat = data.frame(notcoexpressed_mat)
colnames(notcoexpressed_mat) = c("inside","<10kb","<100kb","<1Mb")
notcoexpressed_mat$row = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000")
notcoexpressed_mat$row <- factor(coexpressed_mat$row,levels = c("0-\n10","10-\n20","20-\n50","50-\n100","100-\n200","200-\n500","500-\n1000"))
notcoexpressed_table <- melt(notcoexpressed_mat, id.vars = "row")
colnames(notcoexpressed_table) = c("Gene","CRD","value")


OR_less_10kb = compute_OR(0,1e04)
OR_less_20kb = compute_OR(1e04,2e04)
OR_less_50kb = compute_OR(2e04,5e04)
OR_less_100kb = compute_OR(5e04,1e05)
OR_less_200kb = compute_OR(1e05,2e05)
OR_less_500kb = compute_OR(2e05,5e05)
OR_less_1Mb = compute_OR(5e05,1e06)
OR_all = c(OR_less_10kb[[1]],OR_less_20kb[[1]],OR_less_50kb[[1]],OR_less_100kb[[1]],OR_less_200kb[[1]],OR_less_500kb[[1]],OR_less_1Mb[[1]])
pval_all = c(OR_less_10kb[[2]],OR_less_20kb[[2]],OR_less_50kb[[2]],OR_less_100kb[[2]],OR_less_200kb[[2]],OR_less_500kb[[2]],OR_less_1Mb[[2]])

pdf("Co-expression_analysis_CRD_gene.pdf",paper="a4r")

ggplot(coexpressed_table, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Coexpressed genes") + ylim(0, 65) +
  geom_bar(stat = "identity") +
  geom_text(aes(y=rep(apply(coexpressed_mat[,1:4],1,sum),4),label = c(paste0("(",round(OR_all),")"),rep("",21))),vjust = -.5,  size=6) + labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))

ggplot(notcoexpressed_table, aes(x = Gene, y = value,fill=CRD))+ ggtitle("Not Coexpressed genes") + ylim(0, 10) +
    geom_bar(stat = "identity") +
    labs(x = "Distance (kb)",y = "Fraction associated with the same CRD (%)") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20))

saveRDS(coexpressed_table, file = "coexpressed_table_70.rds")
saveRDS(notcoexpressed_table, file = "notcoexpressed_table_70.rds")

dev.off()

#HiC data, cordata already has Hi-C value computed for the coexpressed pair
compute_ratio_hic <-function(coexpressed=T,CRDmindist,CRDmaxdist,mindist,maxdist,pval.cutoff=0.01,threshold=5){
    if(coexpressed){
      # if the cordata (gene expression) is in the 1Mb distance, and the Hi-C value (cordata$hic) is above the threshold and the pvaluze under the cutoff 
      # and for down: cordata$hic>threshold does not apply
        up   = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$hic>threshold & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
        down = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$pval.adj<pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    } else {
        up = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$hic>threshold & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
        down = sum(cordata$sameCRD>=CRDmindist & cordata$sameCRD<CRDmaxdist & cordata$pval.adj>=pval.cutoff & abs(cordata$distance)>=mindist & abs(cordata$distance)<maxdist,na.rm=T)
    }
    up/down*100
}

#HiC support by gene to gene distance

gene_dist_bins = c(0,1e04,2e04,5e04,1e05,2e05,5e05,1e06,3e09)
coexpressed_mat_hic = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=2)
notcoexpressed_mat_hic = matrix(0,nrow=(length(gene_dist_bins)-1),ncol=2)

for(i in 1:(length(gene_dist_bins)-1)){
        coexpressed_mat_hic[i,1] = compute_ratio_hic(T,0,1e09,gene_dist_bins[i],gene_dist_bins[i+1])
        coexpressed_mat_hic[i,2] = compute_ratio_hic(T,-1,0,gene_dist_bins[i],gene_dist_bins[i+1])
}

# fig 4.4

toplot2 = data.frame(NEU=c(35.95238, 38.22394, 36.69468, 43.73297, 46.57360, 43.96843, 63.79310, 77.17842),
                     MON=c(29.59184, 31.05023, 38.69732, 35.12476, 39.72366, 50.15480, 48.20847, 49.00662),
                     TCL=c(30.39216, 28.00000, 24.82759, 21.48148, 31.85841, 47.48201, 49.41176, 70.68966),
                     dist=c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))


toplot2.molten = melt(toplot2,id.vars="dist")
toplot2.molten$dist = factor(toplot2.molten$dist,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","100-200kb","200-500kb","0.5-1Mb"))
colnames(toplot2.molten)[2] = "CellType"

g <- ggplot(toplot2.molten, aes(x = dist, y = value,fill=CellType))+ ggtitle("HiC support for gene pairs associated with CRD") +
  geom_bar(position="dodge", stat="identity") + theme_classic() + scale_fill_manual(values =  c("#F1BB7B", "#FD6467", "#5B1A18")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+ 
  theme(text = element_text(size=18),axis.title = element_text(size = 15),axis.text = element_text(size = 15),axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Distance between genes",y = "Fraction with HiC support (%)") 
print(g)


pdf("HiC_support_gene_pairs_with_same_CRD.pdf")

toplot = data.frame(counts = coexpressed_mat_hic[,1],dist = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))
toplot$dist <- factor(toplot$dist ,levels = c("0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb",">1Mb"))

ggplot(toplot, aes(x = dist, y = counts))+ ggtitle("HiC support for gene pairs associated with CRD") +
   geom_bar(stat = "identity",fill="#56B4E9") +
   labs(x = "Distance between genes",y = "Fraction with HiC support (%)") +
   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


#HiC support by gene to CRD distance

CRD_dist_bins = c(0,1,1e04,2e04,5e04,1e05,2e05,5e05,1e06)
coexpressed_mat_hic_crddist = matrix(0,nrow=(length(CRD_dist_bins)-1),ncol=1)

for(i in 1:(length(CRD_dist_bins)-1)){
        coexpressed_mat_hic_crddist[i,1] = compute_ratio_hic(T,CRD_dist_bins[i],CRD_dist_bins[i+1],0,1e09)
}

pdf("HiC_support_gene_pairs_with_same_CRD_geneCRD_dist.pdf")
toplot = data.frame(counts = coexpressed_mat_hic_crddist[,1],dist = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb"))
toplot$dist <- factor(toplot$dist ,levels = c("inside","0-10kb","10-20kb","20-50kb","50-100kb","0.1-0.2Mb","0.2-0.5Mb","0.5-1Mb"))

ggplot(toplot, aes(x = dist, y = counts))+ ggtitle("HiC support for gene pairs associated with CRD") +
   geom_bar(stat = "identity",fill="#56B4E9") +
   labs(x = "Distance between genes",y = "Fraction with HiC support (%)") +
   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(text = element_text(size=18),axis.title = element_text(size = 20),axis.text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

save.image('Gene_CRD_associations.RDAta')
