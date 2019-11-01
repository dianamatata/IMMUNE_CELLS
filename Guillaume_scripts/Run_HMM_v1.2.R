library(corrplot)
library(ggplot2)
library(RcppHMM)

chroms = c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22')

datafolder = 'THREE_CELL_TYPES/CLOMICS/'
cell_type_ID1 = 'EGAD00001002670_CLOMICS_v3.0'
cell_type_ID2 = 'EGAD00001002672_CLOMICS_v3.0'
cell_type_ID3 = 'EGAD00001002673_CLOMICS_v3.0'
cell_types = c(cell_type_ID1,cell_type_ID2,cell_type_ID3)
correlation_method = 's'
winsizes = c(1,2,3,5,10)
doHMM = T

for(winsize in winsizes){

##################################
#HMM training and decoding
##################################
if(doHMM){
    for(cell_type_ID in cell_types) {
      segmentation = c(0,1,2,5,10,20,50,100,200,500)
      seglabels = LETTERS[1:(length(segmentation)-1)]

      com1 = paste0('awk \'FNR>1||NR==1\' ',datafolder,cell_type_ID,'/win_',winsize,'/peak_pvalue_win_',winsize,'_chr* > ',datafolder,cell_type_ID,'/win_',winsize,'/peak_pvalue_win_',winsize,'_all.txt')
      cat(com1,'\n')
      system(com1)

      myprofile = read.table(paste0(datafolder,cell_type_ID,'/win_',winsize,'/peak_pvalue_win_',winsize,'_all.txt'),header=T,stringsAsFactors=F)[,2]

      myprofile.discrete = cut_interval(myprofile,n=4,breaks = segmentation,labels = seglabels)
      set.seed(1)
      n <- c("1", "0")
      m <- seglabels
      A <- matrix(c(0.99,0.01,0.05,0.95),nrow = 2,byrow = TRUE)
      B <- matrix(c(1:length(seglabels),length(seglabels):1),nrow = 2,byrow = TRUE)
      B[1,] = B[1,]/sum(B[1,])
      B[2,] = B[2,]/sum(B[2,])
      Pi <- c(0.01,0.99)
      params <- list( "Model" = "HMM",
                          "StateNames" = n,
                          "ObservationNames" = m,
                          "A" = A,
                          "B" = B,
                          "Pi" = Pi)

      HMM <- verifyModel(params)


      HMMFit <- learnEM(HMM,myprofile.discrete,iter=1000,print=F)

      pdf(paste0(datafolder,cell_type_ID,'/win_',winsize,'/HMM_all_chr_win_',winsize,'.pdf'),height=5,width=15)
      for(chr in chroms){
        myprofile = read.table(paste0(datafolder,cell_type_ID,'/win_',winsize,'/peak_pvalue_win_',winsize,'_chr',chr,'.txt'),header=T,stringsAsFactors=F)[,2]
        myprofile.discrete = cut_interval(myprofile,n=4,breaks = segmentation,labels = seglabels)
        predictedStates = as.numeric(forwardBackward(HMMFit,myprofile.discrete))
        write.table(data.frame(x=1:length(predictedStates),state=predictedStates),file=paste0(datafolder,cell_type_ID,'/win_',winsize,'/HMM_CRD_win_',winsize,'_chr',chr,'.txt'),row.names=F,quote=F,sep='\t')
        plot(myprofile,type='l',main=paste0("chr",chr),xlab="Peaks",ylab="-log10(pval)",ylim=c(0,200))
        lines(predictedStates*10+100,col='red')
      }
      dev.off()
    }
}

#############################################################
#Keep only CRDs with at least 2 distinct regulatory elements
#############################################################


for(cell_type_ID in cell_types) {
    mybed = read.table(paste0(datafolder,cell_type_ID,'/merged_residuals.bed.gz'),stringsAsFactors=F)[,c(1:4)]
    mybed = mybed[-which((mybed[,1] == "X" | mybed[,1] == "Y") | mybed[,1] == "M"),]
    mybed[,1] = as.numeric(mybed[,1])
    mybed = mybed[order(mybed[,1]),]
    colnames(mybed) = c("chr","start","end","id")
    nCRDs = 0
    for(chr in rev(chroms)){
        mybed.chr = mybed[which(mybed$chr == as.numeric(chr)),]
        mybed.chr = mybed.chr[order(mybed.chr$start),]
        rownames(mybed.chr) = c(1:nrow(mybed.chr))

        myhmm = read.table(paste0(datafolder,cell_type_ID,'/win_',winsize,'/HMM_CRD_win_',winsize,'_chr',chr,'.txt'),header=T,stringsAsFactors=F)

        CRD_npeaks = rle(myhmm[,2])
        current = 1
        if(myhmm[1,2] == 1){
            firstCRDidx = 1
        } else {
            firstCRDidx = 2
            current = current + CRD_npeaks$lengths[1]
        }

        for(n in seq(firstCRDidx,length(CRD_npeaks$lengths)-1,by=2)){
            # cat("Current peak: ",current,"\n")
            corrIDlist = current:(current+CRD_npeaks$lengths[n]-1)

            got_a_crd = FALSE
            if(length(corrIDlist)>1){
                for(i in 1:(length(corrIDlist)-1)) {
                    for(j in (i+1):length(corrIDlist)){
                        posI = corrIDlist[i]
                        posJ = corrIDlist[j]
                        # cat("i: ",posI,"j: ",posJ,"\n")
                        if(length(intersect(c(mybed.chr[posI,2]:mybed.chr[posI,3]),c(mybed.chr[posJ,2]:mybed.chr[posJ,3]))) == 0){
                            got_a_crd = TRUE
                        }
                    }
                }
            }
            if(!got_a_crd){
                myhmm[current:(current+CRD_npeaks$lengths[n]-1),2] = 0
                # cat("CRD excluded, no peaks:",length(corrIDlist),"\n")
            } else {
                nCRDs = nCRDs+1
                cat("CRD no ",nCRDs,"\r")
            }
            current = current + CRD_npeaks$lengths[n] + CRD_npeaks$lengths[n+1]
        }
        cat("CRD no ",nCRDs,"\n")
        write.table(myhmm,file=paste0(datafolder,cell_type_ID,'/win_',winsize,'/HMM_CRD_FILTERED_win_',winsize,'_chr',chr,'.txt'),row.names=F,quote=F,sep='\t')
    }

    myhmm.all = c()
    myscore = c()
    for(chr in chroms){
        myhmm.all = c(myhmm.all,read.table(paste0(datafolder,cell_type_ID,'/win_',winsize,'/HMM_CRD_FILTERED_win_',winsize,'_chr',chr,'.txt'),header=T,stringsAsFactors=F)[,2],0)
        myscore = c(myscore,read.table(paste0(datafolder,cell_type_ID,'/win_',winsize,'/peak_pvalue_win_',winsize,'_chr',chr,'.txt'),header=T,stringsAsFactors=F)[,2],0)
    }

    pdf(paste0(datafolder,cell_type_ID,'/win_',winsize,'/Peak_length_vs_intensity_FILTERED_win_',winsize,'.pdf'),height=5,width=5)

    CRD_npeaks = rle(myhmm.all)
    CRD_strength = c()
    current = 1
    for(j in seq(2,length(CRD_npeaks$lengths),by=2)){
      current = current + CRD_npeaks$lengths[j-1]
      CRD_strength = c(CRD_strength,mean(myscore[current:(current+CRD_npeaks$lengths[j]-1)]))
      current = current + CRD_npeaks$lengths[j]
    }

    plot(CRD_npeaks$lengths[CRD_npeaks$values==1]+1,CRD_strength,xlab="Number of peaks per CRD",ylab="Average score")
    mycor = as.numeric(cor.test(CRD_npeaks$lengths[CRD_npeaks$values==1],CRD_strength,method='k')$estimate)
    text(max(CRD_npeaks$lengths[CRD_npeaks$values==1])/2,max(CRD_strength)/2,paste("Kendall's corr =",round(mycor,3)))
    hist(CRD_npeaks$lengths[CRD_npeaks$values==1]+1,xlab="Number of peaks per CRD",main=paste0("Total number of CRDs: ",length(CRD_npeaks$lengths[CRD_npeaks$values==1])),breaks=seq(0,150,by=1),xlim=c(0,50))
    dev.off()

    write.table(data.frame(x=1:length(myhmm.all),state=myhmm.all),file=paste0(datafolder,cell_type_ID,'/win_',winsize,'/HMM_CRD_FILTERED_win_',winsize,'_all.txt'),row.names=F,quote=F,sep='\t')

}
}
