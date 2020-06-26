import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from varname import nameof
import numpy as np
import collections
import seaborn as sns
import matplotlib.pyplot as plt


xl = pd.ExcelFile(
    '/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/results/gorilla_panther_grofiler2_pathway_enrichment.xlsx')
print(xl.sheet_names)
gorilla = xl.parse("Gorilla_3diffbackground_3infos")
gorilla = gorilla[gorilla.columns.drop("plots")]
gorilla = gorilla[gorilla.columns.drop("Genes")]
gorilla = gorilla[gorilla['cluster'].isna() != True]
gorilla.shape
go = gorilla[gorilla['FDR q-value'] <= 0.01]
go.shape
go005 = gorilla[gorilla['FDR q-value'] <= 0.05]

subgo = go[['GO term', 'Description', 'FDR q-value', 'type', 'cluster']]
subgo.to_excel(
    '/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/stats.xlsx', index=False)

counting = go.groupby('cluster').count()['Description']
counting_type = go.groupby('type').count()

counting_005 = go005.groupby('cluster').count()['Description']





go70 = go[go['cluster'].str.contains('70')]
go70.cluster.unique()
pvalue_threshold = 0.001

go.query('Gender=="Male" & Year=="2014" ')
groups = go.groupby(['cluster', 'type'])
a=groups.size()

df = go.pivot_table(index='cluster', columns='type', values='fare', aggfunc=np.median)
sns.heatmap(df, annot=True, fmt=".1f")
plt.show()

##### to plot genes and crds per hubs
# sns.jointplot(data=titanic, x='age', y='fare', kind='reg', color='g')
# plt.show()
