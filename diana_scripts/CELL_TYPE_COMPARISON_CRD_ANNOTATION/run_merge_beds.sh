/software/UHTS/Analysis/BEDTools/2.26.0/bin/bedtools merge -i ../EGAD00001002672_CLOMICS_v3.0/EGAD00001002672_ALL.ALLchr.module.bed > EGAD00001002672_ALL.ALLchr.module.merged.bed
/software/UHTS/Analysis/BEDTools/2.26.0/bin/bedtools merge -i ../EGAD00001002673_CLOMICS_v3.0/EGAD00001002673_ALL.ALLchr.module.bed > EGAD00001002673_ALL.ALLchr.module.merged.bed
/software/UHTS/Analysis/BEDTools/2.26.0/bin/bedtools merge -i ../EGAD00001002670_CLOMICS_v3.0/EGAD00001002670_ALL.ALLchr.module.bed > EGAD00001002670_ALL.ALLchr.module.merged.bed

cat EGAD00001002672_ALL.ALLchr.module.merged.bed EGAD00001002673_ALL.ALLchr.module.merged.bed | sort -V -k1,1 -k2,2n > 72_73.ALLchr.module.bed
cat EGAD00001002670_ALL.ALLchr.module.merged.bed EGAD00001002673_ALL.ALLchr.module.merged.bed | sort -V -k1,1 -k2,2n > 70_73.ALLchr.module.bed
cat EGAD00001002670_ALL.ALLchr.module.merged.bed EGAD00001002672_ALL.ALLchr.module.merged.bed | sort -V -k1,1 -k2,2n > 70_72.ALLchr.module.bed
cat EGAD0000100267?_ALL.ALLchr.module.merged.bed | sort -V -k1,1 -k2,2n > 70_72_73.ALLchr.module.bed

/software/UHTS/Analysis/BEDTools/2.26.0/bin/bedtools merge -i 70_73.ALLchr.module.bed > 70_73.ALLchr.module.merged.bed
/software/UHTS/Analysis/BEDTools/2.26.0/bin/bedtools merge -i 70_72.ALLchr.module.bed > 70_72.ALLchr.module.merged.bed
/software/UHTS/Analysis/BEDTools/2.26.0/bin/bedtools merge -i 72_73.ALLchr.module.bed > 72_73.ALLchr.module.merged.bed


