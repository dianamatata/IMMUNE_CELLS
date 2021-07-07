#!/bin/bash
# 4_run_homner_quantify.sh

# cell type list and histone modification list
declare -a Cell_Types=("EGAD00001002673" "EGAD00001002672" "EGAD00001002670")
declare -a Histone_Modifications=("H3K4me1" "H3K27ac")

# folders
DATA_FOLDER=/data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells
THREE_CELL_TYPE_FOLDER=/data/unige/funpopgen/davalos/project/IMMUNE_CELLS
OUTPUT_FOLDER=$THREE_CELL_TYPE_FOLDER/step3_homer_quantify

for cell_type in ${Cell_Types[@]}; do
    for histone in ${Histone_Modifications[@]}; do
        echo "$cell_type""_$histone"

        # check the paths not ready to get launched, where are and how do we decide about samples to exclude? METABAM folder?
        cmd="python3 run_homer_quantif.py $THREE_CELL_TYPE_FOLDER/step1_metabam/$histone"".metafile.bam $DATA_FOLDER/$cell_type/$histone/ $OUTPUT_FOLDER/$histone $THREE_CELL_TYPE_FOLDER/samples_to_exclude.txt"
        echo "$cmd"
        #eval $cmd
    done
done
