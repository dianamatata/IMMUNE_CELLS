#!/bin/bash

# Goal: create .bai index files for bam files

DATA_FOLDER=/home/users/a/avalosma/scratch/Blueprint/ChIP-Seq
for cell_type in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
        for histone in 'H3K4me1' 'H3K27ac' ; do
		for file in $DATA_FOLDER/${cell_type}/$histone/*.bam ; do
			cmd="samtools index -b $file"
			sbatch --partition=mono-shared-EL7 --time=03:00:00 --wrap="$cmd"
		done
	done
done
