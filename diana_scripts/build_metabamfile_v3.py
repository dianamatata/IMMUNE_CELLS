import os
import glob
import subprocess
import sys
import random

datafolder = sys.argv[1]
outputfile = sys.argv[2]
outputtag = outputfile.replace('.bam', '')

# random selection of samples indices, making sure no sample is picked twice
number_of_samples = 50
files = glob.glob(datafolder + '/*.bam')  # number_of_files > number_of_samples
samples_selected = list(range(len(files)))
random.shuffle(samples_selected)  # take the first 50 elements of the number of files shuffled

for j in range(number_of_samples):
    i = samples_selected[j]
    # glob.glob cannot find recursively in subfolder but needs the exact folder. otherwise: https://stackoverflow.com/questions/2186525/how-to-use-glob-to-find-files-recursively
    subsample = "/data/unige/funpopgen/odelanea/SHARE/downsampleBAM/bin/sampleBAM " + files[
        i] + " 1000000 " + outputtag + "_" + str(j) + ".bam"
    print(subsample)
    os.system(subsample)

mergebam = '/software/UHTS/Analysis/samtools/1.4/bin/samtools merge ' + outputfile + ' ' + outputtag + "_*.bam"
print(mergebam)
os.system(mergebam)
