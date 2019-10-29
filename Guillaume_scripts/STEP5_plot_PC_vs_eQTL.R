library(stringr)
args = commandArgs(trailingOnly=TRUE)


datafiles = args[1]
# datafiles = "QTLtools_EGAD00001002670_H3K4me1_1000perm/permutations_cov*.significant.txt"
datafolder = str_extract(datafiles, "QTL[A-Za-z0-9_]+")
chosen_PC = as.numeric(args[2])
cat(chosen_PC,args[1],sep='\n')
#chosen_PC = 22
# chosen_PC = chosen_PC

myfiles = Sys.glob(datafiles)
mydat = c()
for(myfile in myfiles){
	cat(myfile,"\n")
	con = pipe(paste("wc -l", myfile))
	mydat = rbind(mydat,scan(con,what='',quiet=T))
	close(con)
}


covariates = as.numeric(sub("cov","",str_extract(mydat[,2], "cov([0-9]+)")))

mydat = mydat[order(covariates),]

covariates = covariates[order(covariates)]

eqtl = mydat[,1]

pdf(paste(datafolder,"/cov_vs_eqtl.pdf",sep=""))
plot(covariates, eqtl,type='b',col='darkblue',xlab="PC",ylab="# of eQTLs",cex.axis=1.5,cex.lab=1.5,bty="n")
abline(v=chosen_PC,lty=2,lwd=2,col='red')
dev.off()
