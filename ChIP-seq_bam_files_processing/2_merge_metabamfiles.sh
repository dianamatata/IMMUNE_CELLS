#!/bin/bash

### Step 2: Once that we have the subsample of histone marks / cell type
# We merge the results of all the cell types, to have a common database of histone marks per histone mark type
# Itr is essential to sort the results in order to call Homer afterwards  

FOLDER1="/home/users/a/avalosma/scratch/step1_metabam"
mkdir -o $FOLDER1

# merge metabamfiles with samtools
samtools merge $FOLDER1/H3K4me1.metafile.bam $FOLDER1/EGAD0000100267?_H3K4me1_50samples.metafile.bam
echo "H3K4me1.metafile.bam created"
samtools merge $FOLDER1/H3K27ac.metafile.bam $FOLDER1/EGAD0000100267?_H3K27ac_50samples.metafile.bam
echo "H3K27ac.metafile.bam created"

# sort metabam files
samtools sort -o $FOLDER1/H3K27ac.sorted.metafile.bam $FOLDER1/H3K27ac.metafile.bam      
samtools sort -o $FOLDER1/H3K4me1.sorted.metafile.bam $FOLDER1/H3K4me1.metafile.bam  
echo "H3K27ac.sorted.metafile.bam and H3K4me1.sorted.metafile.bam created: H3K4me1.metafile.bam and H3K27ac.metafile.bam sorted"
