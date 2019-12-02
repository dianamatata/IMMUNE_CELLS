#!/bin/bash

HOMERPATH=/software/UHTS/Analysis/homer/4.9/bin
METABAMFOLDER=/data/unige/funpopgen/davalos/project/IMMUNE_CELLS/step1_metabam
declare -a Histone_Modifications=("H3K4me1" "H3K27ac")

OUTPUT_FOLDER=/data/unige/funpopgen/davalos/project/IMMUNE_CELLS/step2_peaks
mkdir -p $OUTPUT_FOLDER

metabamtag='ls *.metafile.bam'
echo $metabamtag

for histone in ${Histone_Modifications[@]}; do

    mkdir -p $OUTPUT_FOLDER/$histone
    mkdir -p $OUTPUT_FOLDER/$histone/META

 
    # Make TagDir for metabamfile
    cmd="$HOMERPATH/makeTagDirectory $OUTPUT_FOLDER/$histone/META $METABAMFOLDER/$histone"".metafile.bam"
    echo "$cmd"


    # Meta-sample peak calling for histone marks:
    # If "-o auto" is specified, the peaks will be written to: "<tag directory>/regions.txt" (-style histone)
    cmd="$HOMERPATH/findPeaks $OUTPUT_FOLDER/$histone/META -style histone -o auto"
    echo "$cmd"

    
    # Convert peaks from homer-format to BED format
    for tag in $metabamtag; do    
        com3="$HOMERPATH/pos2bed.pl $OUTPUT_FOLDER/$histone/META/regions.txt > "+$OUTPUT_FOLDER/$histone/META/""+metabamtag+".bed"
        echo "$cmd"
    done


done

