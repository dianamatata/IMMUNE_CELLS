#!/bin/bash

declare -a Histone_Modifications=("H3K4me1" "H3K27ac")


for histone in ${Histone_Modifications[@]}; do
    FILENAME="samtools_merging_and_sorting_$histone"

    # merge metabamfiles with samtools
    cmd="samtools merge ${histone}.metafile.bam EGAD0000100267?_${histone}_50samples.metafile.bam && samtools sort -o ${histone}.sorted.metafile.bam ${histone}.metafile.bam"
    
    job="bsub -o ${FILENAME}.out -R \"rusage[mem=2000]\" -M 2000000 -n 30 -R \"span[hosts=1]\" -J \"$FILENAME\" \"$cmd\""
    
    echo "$job"
    eval $job
done



