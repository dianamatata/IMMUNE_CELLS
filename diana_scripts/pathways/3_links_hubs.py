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


for cell in ('70', '72', '73'):
    CRD_genes_file = '/Users/dianaavalos/Programming/pathways/mapping_gene_CRD_mean_ALL_%s.txt' % (cell)
    CRD_genes_list = pd.read_csv(CRD_genes_file, header=None, sep=' ')
    CRD_cluster_file = '/Users/dianaavalos/Programming/pathways/trans_hubs_crds/dict_%s_*' % (cell)
    crd_clusters_list = glob.glob(CRD_cluster_file, recursive=True)
    CRD_genes_list.columns = ["phenotype_ID", "phenotype_ID_chr", "phenotype_ID_start", "phenotype_ID_end",
                              "phenotype_ID_strand", "nb_variants", "distance", "CRD_ID", "CRD_ID_chr", "CRD_ID_start",
                              "CRD_ID_end", "nominal_pval", "slope", "top_variant"]
    # for Homer
    CRD_genes_list['phenotype_ID_chr'] = CRD_genes_list['phenotype_ID_chr'].apply(lambda x: 'chr' + str(x))

    # loop on clusters
    for cluster_i_file in crd_clusters_list:
        #        print(cluster_i_file)
        crds_in_cluster = pd.read_table(cluster_i_file, sep=',', header=None, usecols=[0])[0]
        genes = pd.DataFrame(columns=["phenotype_ID", "phenotype_ID_chr", "phenotype_ID_start", "phenotype_ID_end",
                                      "phenotype_ID_strand", "CRD_ID"])
        for CRD in crds_in_cluster:
            # CRD = '11_internal_10924'
            subdf = CRD_genes_list.loc[CRD_genes_list['CRD_ID'] == CRD].iloc[:, [0, 1, 2, 3, 4, 7]]
            genes = gene.append(subdf)

        if (len(genes) > 0):
            outname = "/Users/dianaavalos/Programming/pathways/trans_hubs_genes/dict_%s_%s.csv" % (
                cell, cluster_i_file.split("_")[4].split('.')[0])
            genes.to_csv(outname, index=False, sep='\t')

