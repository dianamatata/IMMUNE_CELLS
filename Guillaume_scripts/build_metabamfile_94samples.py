import os
import glob
import subprocess
import sys

datafolder = sys.argv[1]
outputfile = sys.argv[2]
samplelist = sys.argv[3]
outputtag = outputfile.replace('.bam','')
#print(samfile)
#print(datafolder)
#print(outputfile)

i = 1
with open(samplelist, 'r') as input:
	for filename in input:
		filename = filename.rstrip("\n\r")
		subsample = "/data/unige/funpopgen/odelanea/SHARE/downsampleBAM/bin/sampleBAM "+datafolder+filename+" 1000000 "+outputtag+"_"+str(i)+".bam"
		print(subsample)
		os.system(subsample)
		i = i+1

mergebam = '/software/UHTS/Analysis/samtools/1.4/bin/samtools merge '+outputfile+' '+outputtag+"_*.bam"
print(mergebam)
os.system(mergebam)


