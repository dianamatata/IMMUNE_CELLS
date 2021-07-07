import os
import sys
import pandas as pd

mergebed = sys.argv[1]
threshold = sys.argv[2]
outfolder = sys.argv[3]


if not os.path.exists(outfolder):
	os.makedirs(outfolder)
if not os.path.exists(outfolder+'/OUT/'):
	os.makedirs(outfolder+'/OUT/')
if not os.path.exists(outfolder+'/ERR/'):
	os.makedirs(outfolder+'/ERR/')

chroms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']

clomicspath = '/data/unige/funpopgen/odelanea/SHARE/clomics/bin/clomics'

print("Starting building trees and calling CRDs")

#run clomics correlations

for chrom in chroms:
	com1 = '\' '+clomicspath+' build --bed '+mergebed+' --out '+outfolder+'/build.'+chrom+' --region '+chrom+' \''
	bsubcom1 = 'bsub -K -o '+outfolder+'/OUT/build_'+chrom+'.out -e '+outfolder+'/ERR/build_'+chrom+'.err -J build_'+chrom+' -q normal -M 16000000 '
	print(bsubcom1+com1)
	os.system(bsubcom1+com1)
	
	com2 = '\' '+clomicspath+' call --tree '+outfolder+'/build.'+chrom+' --threshold '+threshold+' --out '+outfolder+'/call_threshold'+threshold+'.'+chrom+' \''
	bsubcom2 = 'bsub -K -o '+outfolder+'/OUT/call_'+chrom+'.out -e '+outfolder+'/ERR/call_'+chrom+'.err -J call_'+chrom+' -q normal -M 16000000 '
	print(bsubcom2+com2)
	os.system(bsubcom2+com2)
        
	com3 = '\'  '+clomicspath+' topo --bed '+mergebed+' --chr '+chrom+' --tree '+outfolder+'/call_threshold'+threshold+'.'+chrom+' --out '+outfolder+'/topo_threshold'+threshold+'.'+chrom+' \''
	bsubcom3 = 'bsub -o '+outfolder+'/OUT/topo_'+chrom+'.out -e '+outfolder+'/ERR/topo_'+chrom+'.err -J topo_'+chrom+' -q normal -M 16000000 '
	print(bsubcom3+com3)
	os.system(bsubcom3+com3)

