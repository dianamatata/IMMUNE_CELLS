library(reshape2)
library(gplots)
library(colorspace)
library(metap)
library(data.table)

# /srv/nasac.unige.ch/funpopgen/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0

z_statistic<-function(r1,r2,n1,n2){
	r1p = 0.5*log((1+r1)/(1-r1))
	r2p = 0.5*log((1+r2)/(1-r2))
	z = (r1p-r2p)/sqrt(1/(n1-3)+1/(n2-3))
	2*pnorm(-abs(z))
}

t_statistic<-function(r,n){
	t = r/sqrt((1-r^2)/(n-2))
	2*pt(-abs(t),n-2)
}

pval_sumlog<-function(x){
	y = -log10(sumlog(x)$p+1e-300)
	y
}

draw_CRD<-function(CRD_peaks_track, myrange){
	current = 0
	j=1
	n = length(CRD_peaks_track)
	while(j <length(CRD_peaks_track)){
		current = CRD_peaks_track[j]
		if(current == 1){
			CRDstart = j-3
			CRDend = which.min(CRD_peaks_track[j:length(CRD_peaks_track)])+j-3
			segments(CRDstart,0,CRDend,0,lwd=2)
			segments(CRDstart,0, (CRDend+CRDstart)/(2),(CRDend-CRDstart),lwd=2)
			segments(CRDend,0, (CRDend+CRDstart)/(2),(CRDend-CRDstart),lwd=2)
			j = CRDend+1
		}
		j=j+1
	}
}


getcolor <- function(x,n=100){
    sequential_hcl(n, h = 260, c = c(0, 100), l = c(100, 30), power = 6,alpha=0.5)[cut(x,n)]
}



args = commandArgs(trailingOnly=TRUE)

datafolder = args[1]
# datafolder = 'THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_CORR_v2.0/'
n_samples = as.numeric(args[2])
# n_samples = 165
chroms = args[3]
# chroms = '21'

#chroms = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')
#chroms = '22'

winsize = 3
winsize_boundaries = 3
plotrange = 100

for(chr in rev(chroms)){
	cat(paste("Running chromosome",chr),'\n')
  	#mydat.cor = read.table(paste(datafolder,'/corr.',chr,sep=''),sep=' ')[,c(1,2,5,6,7)]
  # added for analysis
  chr='21'
  mydat.cor = read.table('/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/methylomics/corr.21.gz',sep=' ')[,c(1,2,5,6,7)]
  
	colnames(mydat.cor) = c("ID1","ID2","PeakName1","PeakName2","cor")
	mydat.peakID = unique(mydat.cor[,c(1,3)])
	mydat.peakID = mydat.peakID[order(mydat.peakID[,1]),]
	tmp <- mydat.cor[(mydat.cor$ID2-mydat.cor$ID1) <= plotrange,]
	newx = (tmp$ID1+tmp$ID2)/2
	newy = tmp$ID2 - tmp$ID1
	cor <- tmp$cor
	Mtoplot <- data.frame("x"=newx,"y"=newy,"corr"=cor^2)
	nmax = max(mydat.cor$ID1)
	com1 = paste("cat ",datafolder,'/topo_threshold2.',chr," | awk -v OFS='\t' '{if($17==-1) print ",chr,",$18,$19,$4,$1}' > ",datafolder,'/topo.',chr,".bed",sep="")
	cat(com1,"\n")
        system(com1)
	com2 = paste("cat ",datafolder,'/topo_threshold2.',chr," | awk '{if($17==1 && $26>1) print $29,$30}' > ",datafolder,'/CRD_peak_IDs.',chr,".txt",sep="")
        cat(com2,"\n")
	system(com2)
	CRD_peak_table = read.table(paste(datafolder,'/CRD_peak_IDs.',chr,".txt",sep=""),stringsAsFactors=F)
	CRD_peak_IDs = unique(sort(unlist(apply(CRD_peak_table,1,function(x){seq(x[1],x[2],by=1)}))))
	CRD_peaks = read.table(paste(datafolder,'/topo.',chr,'.bed',sep=''),stringsAsFactors=F)
	CRD_peaks_track = vector(mode = "numeric", length = nmax)
	CRD_peaks_track[intersect(CRD_peaks[,5],CRD_peak_IDs)]=1

	mydat.cor = data.table(mydat.cor)
	right_corr = mydat.cor[,.(corr_in_win = list(cor[abs(ID2-ID1)<=winsize])),by = ID1]
	left_corr =  mydat.cor[,.(corr_in_win = list(cor[abs(ID2-ID1)<=(winsize+1)][-1])),by = ID2]
	all_corr = mapply(c,right_corr[,2][[1]],left_corr[,2][[1]])
	all_pvals = sapply(all_corr,t_statistic,n=n_samples)
	enrichment_sym = sapply(all_pvals,pval_sumlog)

	cat("Done!\n")

	write.table(file=paste0(datafolder,"/peak_pvalue_chr",chr,".txt"),data.frame(x=1:nmax,p=enrichment_sym),quote=F,row.names=F)
	# write.table(file=paste0(datafolder,"/boundary_pvalue_chr",chr,".txt"),data.frame(x=1:nmax,p=boundaries),quote=F,row.names=F)

	png(paste(datafolder,'/Correlation_plot_win_',winsize,'_chr',chr,'_v2.3.png',sep=''),height=300,width=2000)
	plot(Mtoplot$x,Mtoplot$y,col=getcolor(Mtoplot$cor),pch=15,xlab='',ylab='',main = paste("Correlation chr",chr,datafolder))
	dev.off()

	png(paste(datafolder,'/Correlation_plot_win_',winsize,'_with_CRD_chr',chr,'_v2.3.png',sep=''),height=300,width=2000)
	plot(Mtoplot$x,Mtoplot$y,col=getcolor(Mtoplot$cor),pch=15,xlab='',ylab='',main = paste("Correlation chr",chr,datafolder))
	draw_CRD(CRD_peaks_track,plotrange)
	dev.off()

	png(paste(datafolder,'/Correlation_plot_win_',winsize,'_with_peaks_chr',chr,'_v2.3.png',sep=''),height=300,width=2000)
	plot(Mtoplot$x,Mtoplot$y,col=getcolor(Mtoplot$cor),pch=15,xlab='',ylab='',main = paste("Correlation chr",chr,datafolder))
	lines(enrichment_sym*0.1+plotrange/2,xlab='Peaks',ylab='-log10(pval)',lwd=2)
	dev.off()

	pdf(paste(datafolder,'/Peak_plot_win_',winsize,'_chr',chr,'_v2.3.pdf',sep=''),height=3,width=10)
	# mycolors = sequential_hcl(100, h = 260, c = c(0, 100), l = c(100, 30), power = 4)
	# heatmap.2(matrix(Mtoplot$cor),col= mycolors,main='',dendrogram='none',col, trace='none',Rowv=F,Colv=F,density.info='density')
	plot(enrichment_sym,xlab='Peaks',ylab='-log10(pval)',main='Peaks',type='l',cex.axis=1.3,cex.lab=1.3)
	# plot(boundaries,xlab='Peaks',ylab='-log10(pval)',main='Boundaries',type='l')
	dev.off()
}
