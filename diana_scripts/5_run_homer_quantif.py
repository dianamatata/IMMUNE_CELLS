import os
import glob
import subprocess
import sys
import re

homerpath = "/software/UHTS/Analysis/homer/4.9/bin/"

metabamfile = sys.argv[1]
datafolder = sys.argv[2]
homerfolder = sys.argv[3]
excluded_samples = sys.argv[4]

outfile_folder = homerfolder+'/OUT'
errfile_folder = homerfolder+'/ERR'

if not os.path.exists(homerfolder):
        os.makedirs(homerfolder)

if not os.path.exists(outfile_folder):
	os.makedirs(outfile_folder)

if not os.path.exists(errfile_folder):
	os.makedirs(errfile_folder)

m = re.search('([A-Za-z0-9_.]+).metafile.bam', metabamfile)
metabamtag = m.group(1) 
print(metabamtag)

m2 = re.search('(EGAD[A-Za-z0-9_.]+)', datafolder)
EGA_ID = m2.group(1)

print(EGA_ID)

f = open(excluded_samples, 'r')
excluded = f.readlines()
excluded = [item.strip() for item in excluded]
print(excluded)
f.close()

#make TagDirs for each sample
i =1
for filename in glob.iglob(datafolder+'/*.bam'):
	toexclude = False
	#print(filename)
	m = re.search('(/[A-Za-z0-9_.]+).bam', filename)
	filetag = m.group(1)
	#print(filetag)
	for x in excluded:
		if re.search(x, filename):
			toexclude = True
	if(not toexclude):
		homer_subfolder = homerfolder+'/'+EGA_ID+filetag
		com4 = "\'"+homerpath+"makeTagDirectory "+homer_subfolder+" "+filename+"\'"
		bsubcom1 = 'bsub -o '+outfile_folder+'/makeTagDir_'+EGA_ID+'_'+metabamtag+'_'+str(i)+'.out -e '+errfile_folder+'/makeTagDir_'+EGA_ID+'_'+metabamtag+'_'+str(i)+'.err -J makeTagDir_'+EGA_ID+'_'+metabamtag+'_'+str(i)+' -q normal '
		com4 =  bsubcom1 + com4
		print(com4)
		os.system(com4)
		i=i+1

# Quantify peaks across all individuals at meta-peak calls:
allsubfolder = homerfolder+"/"+EGA_ID+"/*/"
bsubcom2 = 'bsub -o '+outfile_folder+'/quant_'+metabamtag+'.out -e '+errfile_folder+'/quant_'+metabamtag+'.err -J quant_'+metabamtag+' -q normal -w \'ended(\'makeTagDir_'+EGA_ID+'_'+metabamtag+'_*\')\' -M 32000000 '
com5 = "\'"+homerpath+"annotatePeaks.pl "+homerfolder+"/META/"+metabamtag+".bed hg19 -noann -nogene -size given -d "+allsubfolder+" > "+homerfolder+"/"+EGA_ID+"_"+metabamtag+".quantif.txt\'"
com5 = bsubcom2 +com5
print(com5)
os.system(com5)


