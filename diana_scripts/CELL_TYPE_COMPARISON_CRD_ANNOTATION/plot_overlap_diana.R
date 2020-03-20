library(UpSetR) # Creates visualizations of intersecting sets using a novel matrix design
library(data.table) # like pandas in python
library(corrplot) # graphical display of a correlation matrix, confidence interval


# neut 70, mono 72, tcell 73

#### Method 1: Look at the overlap in the surface covered by the CRDs

# from EGAD0000100267*.ALLchr.module.merged.bed
neut=4839
mono=6271
tcell=3086
# not supposed to be used but has been computed, like 70_73.ALLchr.module.merged.bed
neut_mono_merged=6157
neut_tcell_merged=5389
mono_tcell_merged=6391
all_merged=5995
# with bedtools intersect, like 73_72.ALLchr.module.intersect.bed, in ~/scratch/3_CRD, ran by 10_intersect_module_bed.sh
## problem, not the same values than Guillaume. these are the number of fragments, lets compute the number of bp

neut_mono=4953
neut_tcell=2536
mono_tcell=2966
all= 2254 


expressionInput <- c(neutrophil = (neut-all-(neut_mono-all)-(neut_tcell-all)), 
                     monocyte = (mono-all-(neut_mono-all)-(mono_tcell-all)), 
                     tcell = (tcell-all-(neut_tcell-all)-(mono_tcell-all)), 
                     `neutrophil&monocyte` = neut_mono-all,
                     `neutrophil&tcell` = neut_tcell-all,
                     `monocyte&tcell` = mono_tcell-all, 
                     `neutrophil&monocyte&tcell` = all)

# plot
pdf("Overlap_CRD_beds.pdf",paper='a4r')
upset(fromExpression(expressionInput), order.by = "degree",text.scale=1.8)
upset(fromExpression(expressionInput), order.by = "freq",text.scale=1.8,main.bar.color="orange",matrix.color="blue")
dev.off()

#### Method 2: Look at the shared peaks belonging to CRDs among cell types

peakset_neut = scan('peak_IDs_70.txt',what="")
peakset_mono = scan('peak_IDs_72.txt',what="")
peakset_tcel = scan('peak_IDs_73.txt',what="")

neut = length(peakset_neut)
mono = length(peakset_mono)
tcel = length(peakset_tcel)
all = length(intersect(intersect(peakset_neut,peakset_mono),peakset_tcel))
neut_and_mono = length(intersect(peakset_neut,peakset_mono))
neut_and_tcel = length(intersect(peakset_neut,peakset_tcel))
tcel_and_mono = length(intersect(peakset_tcel,peakset_mono))

expressionInput2 <- c(neutrophil = (neut-all-(neut_and_mono-all)-(neut_and_tcel-all)), monocyte = (mono-all-(tcel_and_mono-all)-(neut_and_mono-all)), tcell = (tcel-all-(neut_and_tcel-all)-(tcel_and_mono-all)), `neutrophil&monocyte` = neut_and_mono-all,
                     `neutrophil&tcell` = neut_and_tcel-all,`monocyte&tcell` = tcel_and_mono-all, `neutrophil&monocyte&tcell` = all)

# plot
pdf("Overlap_CRD_beds2.pdf",paper='a4r')
upset(fromExpression(expressionInput2), order.by = "degree",text.scale=1.8)
upset(fromExpression(expressionInput2), order.by = "freq",text.scale=1.8,main.bar.color="orange",matrix.color="blue")
dev.off()



