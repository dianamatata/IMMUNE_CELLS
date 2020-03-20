#!/bin/bash

DATA_FOLDER=/data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells
SAMPLES_TO_EXCLUDE=/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/samples_to_exclude.txt
FOLDER=/home/users/a/avalosma/scratch/step1_metabam
OUTPUT_FOLDER=/home/users/a/avalosma/scratch/step3_homer_quantify

for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
    for histone in 'H3K4me1' 'H3K27ac' ; do
        echo "$cell_type""_$histone"

        # check the paths not ready to get launched, where are and how do we decide about samples to exclude? METABAM folder?
        cmd="python3 run_homer_quantif.py $FOLDER/step1_metabam/$histone"".metafile.bam $DATA_FOLDER/$cell_type/$histone/ $OUTPUT_FOLDER/$histone $SAMPLES_TO_EXCLUDE"
        echo "$cmd"
        #eval $cmd
    done
done
