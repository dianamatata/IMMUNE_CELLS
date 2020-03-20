#!/bin/bash

# add homer to path in .zshrc after downloading in bin so you only need to call makeTagDirectory

METABAMFOLDER=/data/unige/funpopgen/davalos/project/IMMUNE_CELLS/step1_metabam
OUTPUT_FOLDER=/data/unige/funpopgen/davalos/project/IMMUNE_CELLS/step2_peaks
mkdir -p $OUTPUT_FOLDER

metabamtag='ls *.metafile.bam'
echo $metabamtag

for histone in 'H3K4me1' 'H3K27ac' ; do

    mkdir -p $OUTPUT_FOLDER/$histone
    mkdir -p $OUTPUT_FOLDER/$histone/META
 
    # Make TagDir for metabamfile
    cmd="makeTagDirectory $OUTPUT_FOLDER/$histone/META $METABAMFOLDER/$histone"".metafile.bam"
    echo "$cmd"
    eval $cmd

    # Meta-sample peak calling for histone marks:
    # If "-o auto" is specified, the peaks will be written to: "<tag directory>/regions.txt" (-style histone)
    cmd2="$HOMERPATH/findPeaks $OUTPUT_FOLDER/$histone/META -style histone -o auto"
    echo "$cmd2"
    #eval $cmd2

    # Convert peaks from homer-format to BED format
    for tag in $metabamtag; do    
        com3="$HOMERPATH/pos2bed.pl $OUTPUT_FOLDER/$histone/META/regions.txt > "+$OUTPUT_FOLDER/$histone/META/""+metabamtag+".bed"
        echo "$cmd3"
    done


done

