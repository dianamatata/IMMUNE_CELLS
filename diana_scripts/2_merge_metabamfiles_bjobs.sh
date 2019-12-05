#!/bin/bash

# merge metabamfiles with samtools
cmd="samtools merge H3K4me1.metafile.bam EGAD0000100267?_H3K4me1_50samples.metafile.bam && samtools sort -o H3K4me1.sorted.metafile.bam H3K4me1.metafile.bam"
bsub -o samtools_merging_and_sorting_H3K4me1.out -R "rusage[mem=2000]" -M 2000000 -n 30 -R "span[hosts=1]" -J "samtools_merging_and_sorting_H3K4me1" " $cmd"
echo "H3K4me1.metafile.bam created"


cmd="samtools merge H3K27ac.metafile.bam EGAD0000100267?_H3K27ac_50samples.metafile.bam && samtools sort -o H3K27ac.sorted.metafile.bam H3K27ac.metafile.bam"
bsub -o samtools_merging_and_sorting_H3K27ac.out -R "rusage[mem=200]" -M 200000 -n 30 -R "span[hosts=1]" -J "samtools_merging_and_sorting_H3K27ac" " $cmd"
echo "H3K27ac.metafile.bam created"



