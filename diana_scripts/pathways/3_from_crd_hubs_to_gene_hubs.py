import glob
import pandas as pd
import csv


def write_to_file(dataframe, outputfile):
    """
    write dataframe in cvs file
    :param dataframe:
    :param outputfile: cvs filename
    """
    w = csv.writer(open(outputfile, "w"))
    for key, val in dataframe.items():
        w.writerow([key, val])


path = '/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways'
path2='/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/trans'

# with opmmunity file from igraph community algo

for cell in ('70', '72', '73'):
    CRD_genes_list = pd.read_csv("%s/cisCRD-gene_associations/mapping_gene_CRD_mean_ALL_%s.txt" % (path, cell), header=None, sep=' ')
    CRD_cluster_file = pd.read_csv('%s/TRH_%s.txt' % (path2, cell), header=0, sep='\t')
    CRD_genes_list.columns = ["phenotype_ID", "phenotype_ID_chr", "phenotype_ID_start", "phenotype_ID_end",
                              "phenotype_ID_strand", "nb_variants", "distance", "CRD_ID", "CRD_ID_chr", "CRD_ID_start",
                              "CRD_ID_end", "nominal_pval", "slope", "top_variant"]
    # loop on clusters
    for cluster in range(1, max(CRD_cluster_file['communities.membership'])+1):
        CRDS_in_TRH=CRD_cluster_file['communities.names'][CRD_cluster_file['communities.membership'] == cluster]

        genes = pd.DataFrame(columns=["phenotype_ID", "phenotype_ID_chr", "phenotype_ID_start", "phenotype_ID_end",
                                      "phenotype_ID_strand", "CRD_ID"])
        for CRD in CRDS_in_TRH:
            #print(CRD)  # CRD = '11_internal_10924'
            subdf = CRD_genes_list.loc[CRD_genes_list['CRD_ID'] == CRD].iloc[:, [0, 1, 2, 3, 4, 7]]
            genes = genes.append(subdf)

        if (len(genes) > 0):
            outname = "%s/trans_hubs_genes/genes_%s_cluster%s.csv" % (path, cell, cluster)
            genes.to_csv(outname, index=False, sep='\t')
            genes.phenotype_ID.to_csv("%s/trans_hubs_genes/genesonly_%s_cluster%s.csv" % (path, cell, cluster), index=False, sep='\t')
