import glob
import pandas as pd
import csv
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

path = '/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways'
path2 = '/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/trans'


def plot_bars(crd_counts, gene_counts, labels, cell_name, cell):
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars
    fig, ax = plt.subplots()
    rects1 = ax.bar(x - width / 2, crd_counts, width, label='CRDs_TRH')
    rects2 = ax.bar(x + width / 2, gene_counts, width, label='genes_TRH')
    ax.set_ylabel('Counts')
    ax.set_title('Counts of CRDs and Genes in TRH for %s' % cell_name[cell])
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    def autolabel(rects):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = rect.get_height()
            ax.annotate('{}'.format(height),
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom')

    autolabel(rects1)
    autolabel(rects2)
    fig.tight_layout()
    plt.axis('off')
    plt.show()


stats = pd.DataFrame(columns=["TRH_name", "CRDs_TRH", "genes_TRH"])
stats_cells = {}
cell_name = {"70": "Neutrophils", "72": "Monocytes", "73": "T-cells"}

for cell in ('70', '72', '73'):
    sc = pd.DataFrame(columns=["TRH_name", "CRDs_TRH", "genes_TRH"])
    CRD_TRH_file = pd.read_csv('%s/TRH_%s.txt' % (path2, cell), header=0, sep='\t')
    for cluster in range(1, max(CRD_TRH_file['communities.membership']) + 1):  # loop on TRH
        CRDs_TRH = CRD_TRH_file['communities.names'][CRD_TRH_file['communities.membership'] == cluster]
        genes_TRH_file = pd.read_csv("%s/trans_hubs_genes/genesonly_%s_cluster%s.csv" % (path, cell, cluster),
                                     header=None, sep='\t')
        subdf = pd.DataFrame(
            {'TRH_name': ['%s_%s' % (cell, cluster)], 'CRDs_TRH': [len(CRDs_TRH)], 'genes_TRH': [len(genes_TRH_file)]})
        stats = stats.append(subdf)
        sc = sc.append(subdf)  # not working, use dict
        stats_cells[cell] = sc

# plots
for cell in ('70', '72', '73'):
    subset = stats_cells[cell].iloc[0:10]
    labels = subset.TRH_name.values
    crd_counts = subset.CRDs_TRH.values
    gene_counts = subset.genes_TRH.values
    plot_bars(crd_counts, gene_counts, labels, cell_name, cell)
