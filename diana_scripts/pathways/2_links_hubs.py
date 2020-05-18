import os
import pandas as pd
import sys
import re
import csv

CRED = '\033[91m'
CEND = '\033[0m'

inFolder = '/home/users/a/avalosma/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_crds2'
inFolder = '/Users/AvalosDiana/GitHub/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_crds2'


def find_indices_of_occurence(dataframe, crds):
    """
    find in dataframe the indexes in which the crds selected appear
    :param dataframe:
    :param crds:
    dataframe=links
    crds=crd_list
    """
    list_oc = dataframe.index[dataframe.apply(lambda x: x.astype(str).str.contains(crds[0]).any(), axis=1)].tolist()
    for crd in crds[1:]:
        list_oc.extend(dataframe.index[dataframe.apply(lambda x: x.astype(str).str.contains(crd).any(), axis=1)])
    list_oc = list(dict.fromkeys(list_oc))
    return list_oc


def update_dict_node(dict_node, idx_crd, crd_list):
    # add the crds to the dictionary, update the count of occurrence
    # if it is a new crd, add it to the new_crd_list
    new_crd_list = []
    for idx in idx_crd:
        crds = [links.iloc[idx]['id1'], links.iloc[idx]['id2']]
        for crd in crds:
            if crd in dict_node.keys():
                dict_node[crd] += 1
            else:
                dict_node[crd] = 1
                new_crd_list.append(crd)
    new_crd_list = [x for x in new_crd_list if x not in crd_list]
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


def update_links(links, new_links):
    links = new_links.reset_index()
    links = links.drop(['index'], axis=1)
    return links


for cell in ('70', '72', '73'):
    links_path = "%s/LINKS_%s.txt" % (inFolder, cell)
    print(links_path)
    links = pd.read_csv(links_path, sep=' ')
    links = links.reset_index()
    links = links.drop(['index'], axis=1)
    links = links.head(n=500)  # temporary
    node = 0

    while len(links) > 1:
        index = 0
        node += 1
        dict_node = {}
        crd_list = [links.iloc[index]['id1'], links.iloc[index]['id2']]
        idx_crd = find_indices_of_occurence(links, crd_list)  # look for indices of each appearance of the dict entries

        if len(idx_crd) > 0:

            while len(idx_crd) > 0:  # weird statement
                dict_node, crd_list = update_dict_node(dict_node, idx_crd, crd_list)  # add them to the dict
                print('nbr idx: %d, new crds: %d, length links= %d ' % (len(idx_crd), len(crd_list), len(links)))
                temp_links = links
                new_links = links.drop(links.index[idx_crd])  # pop rows from dataframe
                links = update_links(links, new_links)
                if (len(crd_list) > 0):
                    idx_crd = find_indices_of_occurence(links, crd_list)
                else:
                    idx_crd = []
                    print('  no new crd')

            # add to structure holding all the dictionaries. weird i reinitiqte it all the time
            dict_sorted = {k: v for k, v in sorted(dict_node.items(), key=lambda item: item[1], reverse=True)}
            file = "%s/dict_%s_%d.csv" % (inFolder, cell, node)
            print(CRED + "links=%d" % len(links) + CEND)
            print(CRED + "dict_%s_%d.csv" % (cell, node) + CEND)
            with open(file, 'wb') as f:
                w = csv.writer(f)
                w.writerows(dict_sorted.items())

        else:
            print('   else loop')
            links = links.drop(links.index[0])
            links = new_links.reset_index()
            links = links.drop(['index'], axis=1)
            print('no clusters, len links %d' % len(links))
            idx_crd = find_indices_of_occurence(links, crd_list)
