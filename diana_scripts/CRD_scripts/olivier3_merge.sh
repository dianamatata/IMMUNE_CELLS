#!/bin/bash

for c in $(seq 1 22); do
	LDIST=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/LCL_ALL.chr$c\.distance.txt.gz
	FDIST=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/FIB_ALL.chr$c\.distance.txt.gz
	SDIST=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/SUB_ALL.chr$c\.distance.txt.gz
	paste -d' ' <(zcat $LDIST) <(zcat $SDIST | cut -d' ' -f12-) <(zcat $FDIST | cut -d' ' -f12-)
done | gzip -c > LSF_ALL.chrALL.distance.txt.gz

for c in $(seq 1 22); do
	LSUB=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/LCL_ALL.chr$c\.subset.txt.gz
	FSUB=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/FIB_ALL.chr$c\.subset.txt.gz
	SSUB=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/SUB_ALL.chr$c\.subset.txt.gz
	paste -d' ' <(zcat $LSUB) <(zcat $SSUB | cut -d' ' -f12-) <(zcat $FSUB | cut -d' ' -f12-)
done | gzip -c > LSF_ALL.chrALL.subset.txt.gz

