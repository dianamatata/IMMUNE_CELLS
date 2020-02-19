import os
import glob
import subprocess
import sys
import random
from pathlib import Path

datafolder = sys.argv[1]
outputfile = sys.argv[2]
samples_to_exclude_file = sys.argv[3]
outputtag = outputfile.replace('.bam', '')

## next is for testing
# samples_to_exclude_file= 'samples_to_exclude.txt'
# bam_names='list_bam_EGAD00001002672_H3K27ac.txt'
# bam_files = open(bam_names, "r").readlines()
# bam_files=[x.replace('\n', '') for i,x in enumerate(bam_files)]


"""
instead of picking the first bam_files, random selection while making sure no sample is picked twice 
by getting the index of the bam_files in a list and shuffling it and taking the num_bam_files first indices
and also removing the indices of the file not conform 
"""

samples_to_exclude = open(samples_to_exclude_file, "r").readlines()
samples_to_exclude = [x.replace('\n', '') for i, x in enumerate(samples_to_exclude)]

num_bam_files = 50
bam_files = glob.glob(datafolder + '/*.bam')
bam_files_filtered = [x for i, x in enumerate(bam_files) if
                      Path(x).stem not in samples_to_exclude]  # remove forbidden samples, and also only take the file without full path and suffix '.bam'
random.shuffle(bam_files_filtered)  # randomize order
samples_selected = list(range(len(bam_files_filtered)))

"""
subsample 50  bam_files (1 bam/individual and take 1000000 samples in each bam file)
keep track of the selected one through list bam_files_filtered, and write the list in file bam_files_selected.txt
"""
# to make glob.glob recursive:   recursive=True or  https://stackoverflow.com/questions/2186525/how-to-use-glob-to-find-bam_files-recursively

bam_files_selected = []

for i in range(num_bam_files):
    subsample = "/home/users/a/avalosma/bin/downsampleBAM/bin/sampleBAM " + bam_files_filtered[
        i] + " 1000000 " + outputtag + "_" + str(i+1) + ".bam"
    bam_files_selected.append(bam_files_filtered[i])
    print(subsample)
    os.system(subsample)

with open('bam_files_selected.txt', 'w') as f:  # a is append mode, w rewrite the file
    f.write("\n".join(bam_files_selected))

mergebam = '/software/UHTS/Analysis/samtools/1.4/bin/samtools merge ' + outputfile + ' ' + outputtag + "_*.bam"
print(mergebam)
os.system(mergebam)
