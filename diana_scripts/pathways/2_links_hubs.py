import os
import pandas as pd
import sys
import re
import csv

inFolder='/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_crds2'
inFolder='/Users/AvalosDiana/GitHub/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_crds2'

def find_indices_of_occurence(dataframe, crds):
    """
    find in dataframe the indexes in which the crds selected appear
    :param dataframe:
    :param crds:
    """
    list = dataframe.index[dataframe.apply(lambda x: x.astype(str).str.contains(crds[0]).any(), axis=1)].tolist()
    for crd in crds[1:]:
        list.extend(dataframe.index[dataframe.apply(lambda x: x.astype(str).str.contains(crd).any(), axis=1)])
    return list


def update_dict_node(dict_node, idx_crd):
    # add the crds to the dictionary, update the count of occurrence
    # if it is a new crd, add it to the new_crd_list
    new_crd_list = []
    for idx in idx_crd:
        # print (links.iloc[idx]) #true index and count of lines, Name: 26, idx=19, reindex?
        crds = [links.iloc[idx]['id1'], links.iloc[idx]['id2']]
        for crd in crds:
            if crd in dict_node.keys():
                dict_node[crd] += 1
            else:
                dict_node[crd] = 1
                new_crd_list.append(crd)
    return dict_node, new_crd_list


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
    links_path = "%s/LINKS_%s.txt" % (inFolder,cell)
    print(links_path)
    links = pd.read_csv(links_path, sep=' ')
    node = 0
    while len(links) > 1:
        # for each node
        for index, row in links.head(n=1).iterrows():
            print(index)
            node += 1
            dict_node = {row['id1']: 1, row['id2']: 1}
            # look for indices of each appearance of the dict entries
            crd_list = [row['id1'], row['id2']]
            while len(crd_list) > 0:
                idx_crd = find_indices_of_occurence(links, crd_list)
                # add them to the dict (+ new crd list) or increase their count
                dict_node, new_crd_list = update_dict_node(dict_node, idx_crd)
                # pop rows from dataframe
                new_links = links.drop(links.index[idx_crd])
                # check new crds, loop till no more new crds
                crd_list = new_crd_list
                links = new_links.reset_index()
                links = links.drop(['index'], axis=1)
            # add to structure holding all the dictionaries
            dict_sorted = {k: v for k, v in sorted(dict_node.items(), key=lambda item: item[1], reverse=True)}
            print(len(dict_sorted))
            dict_sorted.to_csv( "%s/dict_%s_%d.csv" % (inFolder,cell, node),
                                sep="\t")
