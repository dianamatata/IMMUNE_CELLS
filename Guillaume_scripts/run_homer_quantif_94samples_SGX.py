import os
import glob
import subprocess
import sys
import re
import subprocess

homerpath = "/software/UHTS/Analysis/homer/4.9/bin/"

metabamfile = sys.argv[1]
datafolder = sys.argv[2]
homerfolder = sys.argv[3]
included_samples = sys.argv[4]

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

EGA_ID = 'LCL'

print(EGA_ID)

f = open(included_samples, 'r')
included = f.readlines()
included = [item.strip() for item in included]
print(included)
f.close()

#make TagDirs for each sample
i =1
for filename in glob.iglob(datafolder+'/*.bam'):
	toinclude = False
	#print(filename)
	m = re.search('(/[A-Za-z0-9_.]+).bam', filename)
	filetag = m.group(1)
	#print(filetag)
	for x in included:
		if re.search(x, filename):
			toinclude = True
	if(toinclude):
		targetfile = subprocess.run(['readlink', '-f',filename], stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
		print(targetfile)
		homer_subfolder = homerfolder+'/'+EGA_ID+filetag
		if not os.path.exists(homer_subfolder):
			 os.makedirs(homer_subfolder)
		com4 = "\'samtools view "+targetfile+" | gzip -c > "+homer_subfolder+"/tmp.sam.gz && "+homerpath+"makeTagDirectory "+homer_subfolder+" "+homer_subfolder+"/tmp.sam.gz -format sam && rm "+homer_subfolder+"/tmp.sam.gz\'"
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


