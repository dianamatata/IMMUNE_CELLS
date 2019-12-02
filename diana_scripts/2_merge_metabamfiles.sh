#!/bin/bash

# merge metabamfiles with samtools
samtools merge H3K4me1.metafile.bam EGAD0000100267?_H3K4me1_50samples.metafile.bam
echo "H3K4me1.metafile.bam created"
samtools merge H3K27ac.metafile.unsorted.bam EGAD0000100267?_H3K27ac_50samples.metafile.bam
samtools sort -T H3K27actmpunsorted -o H3K27ac.metafile.bam H3K27ac.metafile.unsorted.bam
echo "H3K27ac.metafile.bam created"

# sort metabam files
samtools sort -T H3K27ac_temp -o H3K27ac.sorted.metafile.bam H3K27ac.metafile.bam      
samtools sort -T H3K4me1_temp -o H3K4me1.sorted.metafile.bam H3K4me1.metafile.bam  
echo "H3K27ac.sorted.metafile.bam and H3K4me1.sorted.metafile.bam created: H3K4me1.metafile.bam and H3K27ac.metafile.bam sorted"
