#!/bin/bash

query_enhancers=
CRD1="EGAD00001002670_CRDs_info.MOD1.NRE2.txt.gz"
CRD2="EGAD00001002672_CRDs_info.MOD1.NRE2.txt.gz"
CRD3="EGAD00001002673_CRDs_info.MOD1.NRE2.txt.gz"

bedtools intersect -wa -wb \
    -a $query_enhancers \
    -b $CRD1 $CRD2 $CRD3 \
    -names 70 72 73 \
    -sorted \
    > output_intersect.txt


