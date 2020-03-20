#!/bin/bash

DATADIR3=/home/users/a/avalosma/scratch/3_CRD

# step 1
for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
	bedtools merge -i $DATADIR3/${cell_type}_94s.ALLchr.module.bed > $DATADIR3/${cell_type}_94s.ALLchr.module.merged.bed
done

# step 2
cat $DATADIR3/EGAD00001002672_94s.ALLchr.module.merged.bed $DATADIR3/EGAD00001002673_94s.ALLchr.module.merged.bed | sort -V -k1,1 -k2,2n > $DATADIR3/72_73.ALLchr.module.bed
cat $DATADIR3/EGAD00001002670_94s.ALLchr.module.merged.bed $DATADIR3/EGAD00001002673_94s.ALLchr.module.merged.bed | sort -V -k1,1 -k2,2n > $DATADIR3/70_73.ALLchr.module.bed
cat $DATADIR3/EGAD00001002670_94s.ALLchr.module.merged.bed $DATADIR3/EGAD00001002672_94s.ALLchr.module.merged.bed | sort -V -k1,1 -k2,2n > $DATADIR3/70_72.ALLchr.module.bed
cat $DATADIR3/EGAD0000100267?_94s.ALLchr.module.merged.bed | sort -V -k1,1 -k2,2n > $DATADIR3/70_72_73.ALLchr.module.bed

# step 3
bedtools merge -i $DATADIR3/72_73.ALLchr.module.bed > $DATADIR3/72_73.ALLchr.module.merged.bed
bedtools merge -i $DATADIR3/70_73.ALLchr.module.bed > $DATADIR3/70_73.ALLchr.module.merged.bed
bedtools merge -i $DATADIR3/70_72.ALLchr.module.bed > $DATADIR3/70_72.ALLchr.module.merged.bed
bedtools merge -i $DATADIR3/70_72_73.ALLchr.module.bed > $DATADIR3/70_72_73.ALLchr.module.merged.bed

