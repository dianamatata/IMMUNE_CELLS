#!/bin/bash

# inputs

query_enhancers_dist="meth_dists.txt"
query_enhancers_prox="meth_proxs.txt"

CRD1="CRDinfo70s.txt"
CRD2="CRDinfo72s.txt"
CRD3="CRDinfo73s.txt"

CRD1="EGAD00001002670_CRDs_info.MOD1.NRE2.txt.gz"
CRD2="EGAD00001002672_CRDs_info.MOD1.NRE2.txt.gz"
CRD3="EGAD00001002673_CRDs_info.MOD1.NRE2.txt.gz"


# commands

bedtools intersect -wa -wb \
    -a $query_enhancers_dist \
    -b $CRD1 $CRD2 $CRD3 \
    -names 70 72 73 \
    > output_intersect_dist.txt

bedtools intersect -wa -wb -a $query_enhancers_prox -b $CRD1 $CRD2 $CRD3 -names 70 72 73 > output_intersect_prox.txt

# in order to sort file for bedtools
# cat CRDinfo73.txt | sort -V -k1,1 -k2,2n > CRDinfo73s.txt

# https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html


