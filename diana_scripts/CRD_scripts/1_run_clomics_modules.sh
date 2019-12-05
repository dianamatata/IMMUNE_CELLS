#!/bin/bash

# Olivier's code for CRD
CLOMICs=/data/unige/funpopgen/odelanea/SHARE/clomics/bin/clomics
# Input data folder so far because I haven't generated yet
DATAFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0
# My output folder for CRD 2nd part of code project
OUTPUT_FOLDER=/data/unige/funpopgen/davalos/project/IMMUNE_CELLS/CRDs
#Loop on all cell types:
declare -a Cell_Types=("EGAD00001002673" "EGAD00001002672" "EGAD00001002670")

LI=$DATAFOLDER/merged_residuals.bed.gz

for cell_type in ${Cell_Types[@]}; do
    CELL_TYPE="$cell_type""_ALL"

    for c in $(seq 1); do
        LO=$DATAFOLDER/$CELL_TYPE\.chr$c
        echo "$LO"
        echo "$CLOMICs build --bed $LI --region $c --out $LO\.1.tgz --silent && $CLOMICs topo --bed $LI --tree $LO\.1.tgz --chr $c --out $LO\.tree.txt.gz && rm $LO\.*.tgz" | bsub -o $LO\.log -R "rusage[mem=8000]" -M 8000000
    done
done

#  #exit
#
#  for c in $(seq 1 22); do
#          LT=$DATAFOLDER/EGAD00001002670_ALL\.chr$c\.tree.txt.gz
#          LO=$DATAFOLDER/EGAD00001002670_ALL\.chr$c\.module
#          #echo "$CLOMICs call --tree $LT --threshold 2 --out $LO\.txt.gz" | bsub -o $LO\.log -R "rusage[mem=8000]" -M 8000000
#  done
#
#  #exit
#
#  for c in $(seq 1 22); do
#          LO=$DATAFOLDER/EGAD00001002670_ALL\.chr$c\.module
#          zcat $LO\.txt.gz | awk '{ if ($30 == 1 && $25 > 1) print $4 }'
#  done | gzip -c > $DATAFOLDER/EGAD00001002670_ALL.modules.MOD1.NRE2.txt.gz

#done





