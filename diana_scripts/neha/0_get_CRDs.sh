#!/bin/bash

# filter the CRDs
# col 25 N_REG col 30 MOD, MOD=1 if CRD, N_REG nbr of peaks, at least 2 R

OUT_FOLDER=/home/users/a/avalosma/scratch/1_CRD

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
        for c in $(seq 1 22); do
                LO=$OUT_FOLDER/${cell_type}.chr$c
                zcat ${LO}.module.txt.gz | awk '{ if ($30 == 1 && $25 > 1) print $4"\t"'$c'"\t"$17"\t"$18}'
        done | gzip -c > $OUT_FOLDER/$cell_type\_CRDs_info.MOD1.NRE2.txt.gz
done

