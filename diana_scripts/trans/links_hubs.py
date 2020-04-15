import os
import pandas as pd
import sys
import re
import csv


links_path = "/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/LINKS.txt"


def find_indices_of_occurence(dataframe, crds):
    list = dataframe.index[dataframe.apply(lambda x: x.astype(str).str.contains(crds[0]).any(), axis=1)].tolist()
    for crd in crds[1:]:
        list.extend(dataframe.index[dataframe.apply(lambda x: x.astype(str).str.contains(crd).any(), axis=1)])
    return list


def update_dict_node(dict_node, idx_crd):
    # add the crds to the dictionary, update the count of occurrence
    # if it is a new crd, add it to the new_crd_list
    new_crd_list = []
    for idx in idx_crd:
        #print (links.iloc[idx]) #true index and count of lines, Name: 26, idx=19, reindex?
        crds = [links.iloc[idx]['from'], links.iloc[idx]['to']]
        for crd in crds:
            if crd in dict_node.keys():
                dict_node[crd] += 1
            else:
                dict_node[crd] = 1
                new_crd_list.append(crd)
    return dict_node, new_crd_list


def write_to_file(dict_sorted, node):
    w = csv.writer(open("dict_%d.csv" %node, "w"))
    for key, val in dict_sorted.items():
        w.writerow([key, val])


links = pd.read_table(links_path, sep=' ')
all_dicts = []
node=0
while len(links)>1:
    # for each node
    for index, row in links.head(n=1).iterrows():
        node+=1
        dict_node = {row['from']: 1, row['to']: 1}
        # look for indices of each appearance of the dict entries
        crd_list = [row['from'], row['to']]
        while len(crd_list) > 0:
            idx_crd = find_indices_of_occurence(links, crd_list)
            # add them to the dict (+ new crd list) or increase their count
            dict_node, new_crd_list = update_dict_node(dict_node, idx_crd)
            # pop rows from dataframe
            new_links = links.drop(links.index[idx_crd])
            # check new crds, loop till no more new crds
            crd_list = new_crd_list
            links = new_links.reset_index()
            links=links.drop(['index'], axis=1)
        # add to structure holding all the dictionaries
        dict_sorted={k: v for k, v in sorted(dict_node.items(), key=lambda item: item[1], reverse=True)}
        print(len(dict_sorted))
        write_to_file(dict_sorted, node)
        all_dicts.append(dict_sorted)


# save list of dict to file
import pickle
f = open("all_dicts.pkl","wb")
pickle.dump(all_dicts,f)
f.close()


# save dict to file
# https://pythonspot.com/save-a-dictionary-to-a-file/

