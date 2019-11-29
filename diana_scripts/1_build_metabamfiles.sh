#!/bin/bash

# cell type list and histone modification list
declare -a Cell_Types=("EGAD00001002673" "EGAD00001002672" "EGAD00001002670")
declare -a Histone_Modifications=("H3K4me1" "H3K27ac")

# folders
DATA_FOLDER=/data/unige/funpopgen/2backup/grey2/external_data/Blueprint_immune_cells
OUTPUT_FOLDER=/data/unige/funpopgen/davalos/project/IMMUNE_CELLS/step1_metabam
SAMPLES_EXCLUDED=/data/unige/funpopgen/davalos/project/IMMUNE_CELLS/diana_scripts/samples_to_exclude.txt
PY_SCRIPT=/data/unige/funpopgen/davalos/project/IMMUNE_CELLS/diana_scripts/build_metabamfile_v3.py

mkdir -p $OUTPUT_FOLDER


for cell_type in ${Cell_Types[@]}; do
    for histone in ${Histone_Modifications[@]}; do
        echo "$cell_type""_$histone"
        cmd="bsub -o $OUTPUT_FOLDER/buildmetabamfile_$cell_type""_$histone.out -e $OUTPUT_FOLDER/buildmetabamfile_$cell_type""_$histone.err -J buildmetabamfile_$cell_type""_$histone -q normal -M 64000000 \"python3.5 $PY_SCRIPT $DATA_FOLDER/$cell_type/$histone $OUTPUT_FOLDER/$cell_type""_$histone""_50samples.metafile.bam\" $SAMPLES_EXCLUDED"
        echo "$cmd"
        eval $cmd
        #$cmd does not work because somehow it prints the \ in EGAD00001002670_H3K27ac_50samples.metafile.bam\"
    done
done


