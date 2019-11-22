import os
import glob
import subprocess
import sys
import random

datafolder = sys.argv[1]
outputfile = sys.argv[2]
outputtag = outputfile.replace('.bam', '')


"""
instead of picking the first bam_files, random selection while making sure no sample is picked twice 
by getting the index of the bam_files in a list and shuffling it and taking the num_bam_files first indices
and also removing the indices of the file not conform 
"""

num_bam_files = 50
bam_files = glob.glob(datafolder + '/*.bam')
samples_selected = list(range(len(bam_files)))
random.shuffle(samples_selected)

forbidden_index = [i for i,x in enumerate(bam_files) if 'Blueprint' in x]
samples_selected.remove(forbidden_index[0])

"""
subsample 50  bam_files (1 bam/individual and take 1000000 samples in each bam file)
keep track of the selected one through list bam_files_selected, and write the list in file bam_files_selected.txt
"""
# to make glob.glob recursive:   recursive=True or  https://stackoverflow.com/questions/2186525/how-to-use-glob-to-find-bam_files-recursively

bam_files_selected=[]

for j in range(num_bam_files):
    i = samples_selected[j]
    subsample = "/data/unige/funpopgen/odelanea/SHARE/downsampleBAM/bin/sampleBAM " + bam_files[
        i] + " 1000000 " + outputtag + "_" + str(j) + ".bam"
    bam_files_selected.append(bam_files[i])
    print(subsample)
    os.system(subsample)

with open('bam_files_selected.txt','a') as f:  # a is append mode
    f.write("\n".join(bam_files_selected))

mergebam = '/software/UHTS/Analysis/samtools/1.4/bin/samtools merge ' + outputfile + ' ' + outputtag + "_*.bam"
print(mergebam)
os.system(mergebam)
