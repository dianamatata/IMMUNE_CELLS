#!/bin/bash

BIN=/data/unige/funpopgen/odelanea/SGX/V2/binaries/clomics/bin/clomics

LK27=/data/unige/funpopgen/odelanea/SGX/V2/data_chromatin/bed/LCL_H3K27ac.PC50.resid.bed.gz
LME1=/data/unige/funpopgen/odelanea/SGX/V2/data_chromatin/bed/LCL_H3K4me1.PC50.resid.bed.gz
LME3=/data/unige/funpopgen/odelanea/SGX/V2/data_chromatin/bed/LCL_H3K4me3.PC50.resid.bed.gz
FK27=/data/unige/funpopgen/odelanea/SGX/V2/data_chromatin/bed/FIB_H3K27ac.PC10.resid.bed.gz
FME1=/data/unige/funpopgen/odelanea/SGX/V2/data_chromatin/bed/FIB_H3K4me1.PC10.resid.bed.gz
FME3=/data/unige/funpopgen/odelanea/SGX/V2/data_chromatin/bed/FIB_H3K4me3.PC10.resid.bed.gz
SUBS=/data/unige/funpopgen/odelanea/SGX/V2/data_chromatin/bed/78_LCL_samples.txt

for c in $(seq 1 22); do
	HIC_VAL=/data/unige/funpopgen/odelanea/SGX/V2/data_hic/chr$c\/MAPQG0/chr$c\_5kb.RAWobserved.gz
	HIC_KRN=/data/unige/funpopgen/odelanea/SGX/V2/data_hic/chr$c\/MAPQG0/chr$c\_5kb.KRnorm.gz
	LOUT=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/LCL_ALL.chr$c
	POUT=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/PER2_ALL.chr$c
	FOUT=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/FIB_ALL.chr$c
	SOUT=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/SUB_ALL.chr$c

	echo "$BIN hic --bed $LK27 $LME1 $LME3 --hic $HIC_VAL $HIC_KRN 5000 --null 1.0 --region chr$c --out $POUT\.txt.gz --permute 10 --window 1e6" | bsub -R "rusage[mem=7000]" -M 7000000
	#echo "$BIN hic --bed $LK27 $LME1 $LME3 --hic $HIC_VAL $HIC_KRN 5000 --null 1.0 --region chr$c --out $LOUT\.txt.gz" | bsub -R "rusage[mem=5000]" -M 5000000
	#echo "$BIN hic --bed $FK27 $FME1 $FME3 --hic $HIC_VAL $HIC_KRN 5000 --null 1.0 --region chr$c --out $FOUT\.txt.gz" | bsub -R "rusage[mem=5000]" -M 5000000
	#echo "$BIN hic --bed $LK27 $LME1 $LME3 --hic $HIC_VAL $HIC_KRN 5000 --null 1.0 --region chr$c --include-samples $SUBS --out $SOUT\.txt.gz"  | bsub -R "rusage[mem=5000]" -M 5000000
done

