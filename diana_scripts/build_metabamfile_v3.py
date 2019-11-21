import os
import glob
import subprocess
import sys

datafolder = sys.argv[1]
outputfile = sys.argv[2]
outputtag = outputfile.replace('.bam','')
#print(samfile)
#print(datafolder)
#print(outputfile)

i = 1
for filename in glob.iglob(datafolder+'/*.bam'):
	subsample = "/data/unige/funpopgen/odelanea/SHARE/downsampleBAM/bin/sampleBAM "+filename+" 1000000 "+outputtag+"_"+str(i)+".bam"
	print(subsample)
	os.system(subsample)
	i = i+1
	if(i>50):
		break

mergebam = '/software/UHTS/Analysis/samtools/1.4/bin/samtools merge '+outputfile+' '+outputtag+"_*.bam"
print(mergebam)
os.system(mergebam)


