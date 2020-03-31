#!/bin/bash

# Goal: create .bai index files for bam files

DATA_FOLDER=/home/users/a/avalosma/scratch/Blueprint/RNA_seq
for cell_type in 'EGAD00001002671'  'EGAD00001002674' 'EGAD00001002675' ; do
        for file in $DATA_FOLDER/${cell_type}/*.bam ; do
		cmd="samtools index -b $file"
		sbatch -J bai_$file.job --partition=mono-EL7 --time=01:00:00 --wrap="$cmd"
	done
done
