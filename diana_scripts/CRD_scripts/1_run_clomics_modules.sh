#!/bin/bash

BIN=/data/unige/funpopgen/odelanea/SHARE/clomics/bin/clomics

DATAFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0

LI=$DATAFOLDER/merged_residuals.bed.gz

for c in $(seq 1 22); do
        LO=$DATAFOLDER/EGAD00001002670_ALL\.chr$c
        #echo "$BIN build --bed $LI --region $c --out $LO\.1.tgz --silent && $BIN topo --bed $LI --tree $LO\.1.tgz --chr $c --out $LO\.tree.txt.gz && rm $LO\.*.tgz" | bsub -o $LO\.log -R "rusage[mem=8000]" -M 8000000
done

#exit

for c in $(seq 1 22); do
        LT=$DATAFOLDER/EGAD00001002670_ALL\.chr$c\.tree.txt.gz
        LO=$DATAFOLDER/EGAD00001002670_ALL\.chr$c\.module
        #echo "$BIN call --tree $LT --threshold 2 --out $LO\.txt.gz" | bsub -o $LO\.log -R "rusage[mem=8000]" -M 8000000
done

#exit

for c in $(seq 1 22); do
        LO=$DATAFOLDER/EGAD00001002670_ALL\.chr$c\.module
        zcat $LO\.txt.gz | awk '{ if ($30 == 1 && $25 > 1) print $4 }'
done | gzip -c > $DATAFOLDER/EGAD00001002670_ALL.modules.MOD1.NRE2.txt.gz


