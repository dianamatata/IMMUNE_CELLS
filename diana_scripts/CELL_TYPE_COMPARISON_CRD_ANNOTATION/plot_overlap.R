library(UpSetR)
library(data.table)
library(corrplot)

# 1761 70_72_73.ALLchr.module.intersect.bed
#  3717 70_72.ALLchr.module.intersect.bed
#  2092 70_73.ALLchr.module.intersect.bed
# 2586 72_73.ALLchr.module.intersect.bed

#  4839 EGAD00001002670_ALL.ALLchr.module.merged.bed
#  6271 EGAD00001002672_ALL.ALLchr.module.merged.bed
#  3086 EGAD00001002673_ALL.ALLchr.module.merged.bed



expressionInput <- c(neutrophil = (4839-1761-(3717-1761)-(2092-1761)), monocyte = (6271-1761-(3717-1761)-(2586-1761)), tcell = (3086-1761-(2092-1761)-(2586-1761)), `neutrophil&monocyte` = 3717-1761,
 `neutrophil&tcell` = 2092-1761,`monocyte&tcell` = 2586-1761, `neutrophil&monocyte&tcell` = 1761)

pdf("Overlap_CRD_beds.pdf",paper='a4r')
upset(fromExpression(expressionInput), order.by = "degree",text.scale=1.8)
upset(fromExpression(expressionInput), order.by = "freq",text.scale=1.8,main.bar.color="orange",matrix.color="blue")
dev.off()


peakset_neut = scan('../EGAD00001002670_CLOMICS_v3.0/peak_IDs.txt',what="")
peakset_mono = scan('../EGAD00001002672_CLOMICS_v3.0/peak_IDs.txt',what="")
peakset_tcel = scan('../EGAD00001002673_CLOMICS_v3.0/peak_IDs.txt',what="")

neut = length(peakset_neut)
mono = length(peakset_mono)
tcel = length(peakset_tcel)
all = length(intersect(intersect(peakset_neut,peakset_mono),peakset_tcel))
neut_and_mono = length(intersect(peakset_neut,peakset_mono))
neut_and_tcel = length(intersect(peakset_neut,peakset_tcel))
tcel_and_mono = length(intersect(peakset_tcel,peakset_mono))

expressionInput <- c(neutrophil = (neut-all-(neut_and_mono-all)-(neut_and_tcel-all)), monocyte = (mono-all-(tcel_and_mono-all)-(neut_and_mono-all)), tcell = (tcel-all-(neut_and_tcel-all)-(tcel_and_mono-all)), `neutrophil&monocyte` = neut_and_mono-all,
 `neutrophil&tcell` = neut_and_tcel-all,`monocyte&tcell` = tcel_and_mono-all, `neutrophil&monocyte&tcell` = all)

pdf("Overlap_CRD_peaks.pdf",paper='a4r')
upset(fromExpression(expressionInput), order.by = "degree",text.scale=1.8)
upset(fromExpression(expressionInput), order.by = "freq",text.scale=c(1.8,1.8,1.8,1.2,1.8,1.8),main.bar.color="orange",matrix.color="blue")
dev.off()

peakset_neut = as.data.frame(fread('../EGAD00001002670_CLOMICS_v3.0/EGAD00001002670_ALL.ALLchr.peak.txt.gz',header=F))
peakset_mono = as.data.frame(fread('../EGAD00001002672_CLOMICS_v3.0/EGAD00001002672_ALL.ALLchr.peak.txt.gz',header=F))
peakset_tcel = as.data.frame(fread('../EGAD00001002673_CLOMICS_v3.0/EGAD00001002673_ALL.ALLchr.peak.txt.gz',header=F))


compare_CRD <- function(query,reference,threshold=0.5){
    CRD_IDs = unique(query$V2)
    n_replicated = 0
    for(j in 1:length(CRD_IDs)){
        current_CRD = CRD_IDs[j]
        current_peaks = query$V1[which(query$V2 == current_CRD)]
        idx_overlap = which(reference$V1 %in% current_peaks)
        if(length(idx_overlap)>0){
            overlapping_peaks_in_same_CRD = sort(table(reference$V2[idx_overlap]),decreasing=T)[1]
        } else {
            overlapping_peaks_in_same_CRD = 0
        }
        if(overlapping_peaks_in_same_CRD/length(current_peaks)>threshold){
            n_replicated = n_replicated + 1
        }
    }
    list(total= n_replicated,fraction=n_replicated/length(CRD_IDs))
}

mono_vs_neut = compare_CRD(peakset_mono,peakset_neut)
mono_vs_tcel = compare_CRD(peakset_mono,peakset_tcel)

neut_vs_mono = compare_CRD(peakset_neut,peakset_mono)
neut_vs_tcel = compare_CRD(peakset_neut,peakset_tcel)

tcel_vs_mono = compare_CRD(peakset_tcel,peakset_mono)
tcel_vs_neut = compare_CRD(peakset_tcel,peakset_neut)

pdf("CRD_pairwise_comparisons_between_cell_types.pdf")
M = matrix(c(1,neut_vs_mono$fraction,neut_vs_tcel$fraction,mono_vs_neut$fraction,1,mono_vs_tcel$fraction,
    tcel_vs_neut$fraction,tcel_vs_mono$fraction,1),ncol=3,byrow=T)
colnames(M) = c("Neutrophils","Monocytes","T cells")
rownames(M) = c("Neutrophils","Monocytes","T cells")
corrplot(M, method = "number",is.corr=F,col = "black",number.cex=1.5,cl.lim = c(0, 1))
corrplot(M,is.corr=F,cl.lim = c(0, 1),p.mat = M,sig.level=-1,insig = "p-value",number.cex=1.5)
dev.off()


compare_CRD_3way <- function(query,reference1,reference2){
    CRD_IDs = unique(query$V2)
    n_replicated = 0
    for(j in 1:length(CRD_IDs)){
        current_CRD = CRD_IDs[j]
        current_peaks = query$V1[which(query$V2 == current_CRD)]
        overlapping_peaks_ref1 = length(intersect(current_peaks,reference1$V1))
        overlapping_peaks_ref2 = length(intersect(current_peaks,reference2$V1))
        if(overlapping_peaks_ref1/length(current_peaks)>0.5 & overlapping_peaks_ref2/length(current_peaks)>0.5){
            n_replicated = n_replicated + 1
        }
    }
    list(total= n_replicated,fraction=n_replicated/length(CRD_IDs))
}


neut_vs_others = compare_CRD_3way(peakset_neut,peakset_mono,peakset_tcel)
mono_vs_others = compare_CRD_3way(peakset_mono,peakset_neut,peakset_tcel)
tcel_vs_others = compare_CRD_3way(peakset_tcel,peakset_neut,peakset_mono)
