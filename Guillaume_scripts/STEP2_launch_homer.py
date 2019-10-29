import os
import glob
import subprocess
import sys
import re

homerpath = "/software/UHTS/Analysis/homer/4.9/bin/"

metabamfile = sys.argv[1]
datafolder = sys.argv[2]
homerfolder = sys.argv[3]

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

#make TagDir for metabamfile
com1 = homerpath+"makeTagDirectory "+homerfolder+"/META "+metabamfile
print(com1)
#os.system(com1)

#Meta-sample peak calling for histone marks:

com2 = homerpath+"findPeaks "+homerfolder+"/META -style histone -o auto"
print(com2)
#os.system(com2)

# convert peaks from homer-format to BED format

com3 = homerpath+"pos2bed.pl "+homerfolder+"/META/regions.txt > "+homerfolder+"/META/"+metabamtag+".bed"
print(com3)
#os.system(com3)

#sys.exit()
#make TagDirs for each sample
i =1
for filename in glob.iglob(datafolder+'/*.bam'):
	#print(filename)
	m = re.search('([A-Za-z0-9_.]+).bam', filename)
	filetag = m.group(1)
	homer_subfolder = homerfolder+"/ALL_"+filetag
	com4 = "\'"+homerpath+"makeTagDirectory "+homer_subfolder+" "+filename+"\'"
	bsubcom1 = 'bsub -o '+outfile_folder+'/makeTagDir_'+metabamtag+str(i)+'.out -e '+errfile_folder+'/makeTagDir_'+metabamtag+str(i)+'.err -J makeTagDir_'+str(i)+' -q normal '
	com4 =  bsubcom1 + com4
	print(com4)
	#os.system(com4)
	i=i+1

# Quantify peaks across all individuals at meta-peak calls:
allsubfolder = homerfolder+"/ALL_*/"
bsubcom2 = 'bsub -o '+outfile_folder+'/quant_'+metabamtag+'.out -e '+errfile_folder+'/quant_'+metabamtag+'.err -J quant_'+metabamtag+' -q normal -w \'ended(\'makeTagDir_*\')\' '
com5 = "\'"+homerpath+"annotatePeaks.pl "+homerfolder+"/META/"+metabamtag+".bed hg19 -noann -nogene -size given -d "+allsubfolder+" > "+homerfolder+"/"+metabamtag+".quantif.txt\'"
com5 = bsubcom2 +com5
print(com5)
#os.system(com5)


