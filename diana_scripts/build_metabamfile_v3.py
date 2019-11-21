import os
import glob
import subprocess
import sys
import random

datafolder = sys.argv[1]
outputfile = sys.argv[2]
outputtag = outputfile.replace('.bam','')
#print(samfile)
#print(datafolder)
#print(outputfile)

# random selection of samples indices, making sure no sample is picked twice
number_of_samples = 50
files = glob.glob(datafolder+'/*.bam')
samples_selected = list(range(number_of_samples))
random.shuffle(samples_selected)


for i in range(len(files)):
    if i in samples_selected:
        # glob.glob cannot find recursively in subfolder but needs the exact folder. otherwise: https://stackoverflow.com/questions/2186525/how-to-use-glob-to-find-files-recursively
        subsample = "/data/unige/funpopgen/odelanea/SHARE/downsampleBAM/bin/sampleBAM " + files[i] + " 1000000 " + outputtag + "_" + str(i) + ".bam"
        print(subsample)
        os.system(subsample)


mergebam = '/software/UHTS/Analysis/samtools/1.4/bin/samtools merge ' + outputfile + ' ' + outputtag + "_*.bam"
print(mergebam)
os.system(mergebam)
