import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from varname import nameof
import numpy as np
import collections

xl = pd.ExcelFile(
    '/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/results/gorilla_panther_grofiler2_pathway_enrichment.xlsx')
print(xl.sheet_names)
gorilla = xl.parse("gorilla")
panther = xl.parse("panther")
gprofiler2 = xl.parse("gprofiler2")

searchfor = ['immune', 'neutrophil', 'leukocyte', 'T cell']
immune_pathways = pd.DataFrame(columns=['pvalue', 'pathway', 'Cell', 'TRHindex', 'MethodName'])

method = gprofiler2

for methodName, method in zip(('gorilla', 'panther', 'gprofiler2'), (gorilla, panther, gprofiler2)):
    print(methodName)
    method['MethodName'] = methodName
    method['Cell'] = 0
    method['TRHindex'] = 0
    method['cluster'].to_string()
    method2 = method[method.iloc[:, 0].notna()]
    method = method2
    method = method[method.columns.drop(list(method.filter(regex='Unnamed')))]
    # a=method[method['cluster'].isna() == True]
    method['TRHindex'] = method['cluster'].apply(lambda x: x.split('r')[1])

    for cell_type in (70, 72, 73):
        # mask = method['cluster'].str.contains(r'%d' % cell_type, na=True)
        method['Cell'][method['cluster'].str.contains(str(cell_type))] = str(cell_type)

    method['immune_pathway'] = 0
    mask2 = method['pathway'].str.contains('|'.join(searchfor))
    method['immune_pathway'][mask2] = 1
    immune_method = method[['pvalue', 'pathway', 'Cell', 'TRHindex', 'MethodName']][mask2]
    immune_pathways = pd.concat([immune_pathways, immune_method], ignore_index=True)

immune_pathways.to_excel(
    '/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/results/immunepathways.xlsx', index=False)

analysis = pd.DataFrame(columns=['Cell', 'immune', 'neutrophil', 'monocyte', 'leukocyte', 'T cell'])
for cell_type in (70, 72, 73):
    a = immune_pathways[immune_pathways['Cell'] == str(cell_type)]['pathway']
    b = pd.DataFrame([cell_type, a.str.contains(r'immune').sum(), a.str.contains(r'neutrophil').sum(),
                      a.str.contains(r'monocyte').sum(), a.str.contains(r'leukocyte').sum(),
                      a.str.contains(r'T cell').sum()])
    b = pd.DataFrame(np.array(b).reshape(-1, len(b)),
                     columns=['Cell', 'immune', 'neutrophil', 'monocyte', 'leukocyte', 'T cell'])
    analysis = pd.concat((analysis, b), ignore_index=True)

c=list(immune_pathways['pathway'].str.strip().str.split('[\W_]+'))
flattened = [val for sublist in c for val in sublist]
flattened = [word.lower() for word in flattened]
counts = collections.Counter(flattened)
counts.most_common(15)

