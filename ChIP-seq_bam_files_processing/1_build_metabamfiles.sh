#!/bin/bash
 
### step 1 build metabamfiles
# we have 3 cell types and 2 histone marks.
# We take a subsample of the histone marks, 1000 of them in 50 bam files chosen at random

# folders
 
DATA_FOLDER="/home/users/a/avalosma/scratch/Blueprint/ChIP-Seq"
OUTPUT_FOLDER="/home/users/a/avalosma/scratch/step1_metabam"
SAMPLES_EXCLUDED="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/samples_to_exclude.txt"
mkdir -p $OUTPUT_FOLDER $OUTPUT_FOLDER/OUT
 
# scripts called

PY_SCRIPT="/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/1_bis_build_metabamfiles.py"
# takes as input: datafolder outputfile samples_to_exclude_file
 
for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
    for histone in 'H3K4me1' 'H3K27ac' ; do

        FILENAME=${OUTPUT_FOLDER}/${cell_type}_${histone}_buildmetabamfile
 
        cmd="python3 $PY_SCRIPT $DATA_FOLDER/$cell_type/$histone $OUTPUT_FOLDER/${cell_type}_${histone}_50samples.metafile.bam $SAMPLES_EXCLUDED"
 
        # command in cluster
        # sbatch -J $FILENAME.job --partition=mono-EL7 --time=30:00:00 -o $OUTPUT_FOLDER/OUT/$FILENAME.out -e $OUTPUT_FOLDER/OUT/$FILENAME.err --wrap="$cmd"
        sbatch -J $FILENAME.job --partition=mono-EL7 --time=30:00:00 --wrap="$cmd"
    done
done

# info about job time consumption:
# with --time=1:00:00 only 2 bam files subsamples are generated
# so aim at 30:00:00, mem=12 488 480K when we have 5 bamfiles subsampled and merged
