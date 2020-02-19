#!/bin/bash
 
### step 1 build metabamfiles
 
# folders
 
DATA_FOLDER="/home/users/a/avalosma/scratch/Blueprint/ChIP-Seq"
OUTPUT_FOLDER="/home/users/a/avalosma/scratch/step1_metabam"
SAMPLES_EXCLUDED="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/samples_to_exclude.txt"
mkdir -p $OUTPUT_FOLDER
 
# scripts called

PY_SCRIPT="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/1_bis_build_metabamfiles.py"
# takes as input datafolder outputfile samples_to_exclude_file
 
# conditions
 
cell_type="EGAD00001002670"
histone="H3K4me1"

FILENAME=${OUTPUT_FOLDER}/${cell_type}_${histone}_buildmetabamfile
 
cmd="python3 $PY_SCRIPT $DATA_FOLDER/$cell_type/$histone $OUTPUT_FOLDER/${cell_type}_${histone}_50samples.metafile.bam $SAMPLES_EXCLUDED"
 
# command in cluster
sbatch -J $FILENAME.job --partition=mono-shared-EL7 --time=24:00:00 -o $FILENAME.out -e $FILENAME.err --wrap="$cmd"


