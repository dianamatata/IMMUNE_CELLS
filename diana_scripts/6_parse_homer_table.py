import pandas as pd
import re
import sys
import os

filename = sys.argv[1]
assaytag = sys.argv[2]
output = filename.replace('.txt','.bed')
phenotype_to_exclude = filename.replace('.txt','.pheno_to_exclude.txt')

mydat = pd.read_table(filename)

#update column names
new_columns = mydat.columns.values
new_columns[0] = 'ID'
regex=re.compile('(S0[A-Z0-9]{4})')
new_columns[7:mydat.shape[1]] = [m.group(1) for l in list(mydat.ix[:,7:mydat.shape[1]]) for m in [regex.search(l)] if m]
mydat.columns = new_columns

#update chr names
chrIDs = list(mydat.ix[:,1])
chrIDs = [w.replace('chr', '') for w in chrIDs]
mydat.ix[:,1]  = chrIDs

mydatbed = mydat.ix[:,1:4].join(mydat.ix[:,[0,0,4]])
mydatbed = mydatbed.join(mydat.ix[:,7:mydat.shape[1]])

new_columns = mydatbed.columns.values
new_columns[0] = '#Chr'
new_columns[3] = 'PhenoID'
new_columns[4] = 'Group_PhenoID'
mydatbed.columns  = new_columns

mydatbed['PhenoID'] = assaytag + mydatbed['PhenoID'].astype(str)
mydatbed['Group_PhenoID'] = assaytag + mydatbed['Group_PhenoID'].astype(str)

#mydatbed = mydatbed.sort_values(by=['#Chr', 'Start'])
mydatbed.to_csv(output,sep='\t',float_format='%.2f',index=False)

os.system('(head -1 '+output+' && tail -n +2 '+output+'| sort -V -k1,1 -k2,2n) | /software/UHTS/Analysis/samtools/1.4/bin/bgzip > '+output+'.gz')
os.system('/software/UHTS/Analysis/samtools/1.4/bin/tabix -p bed '+output+'.gz')

temp = mydatbed.ix[:,7:mydatbed.shape[1]].abs().sum(axis=1) == 0
mydatpheno = mydatbed.ix[temp,3]
mydatpheno.to_csv(phenotype_to_exclude,sep='\t',float_format='%.2f',index=False)

