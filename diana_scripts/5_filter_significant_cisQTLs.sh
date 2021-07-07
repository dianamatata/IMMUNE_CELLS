#!/bin/bash

# after running create_permutations.sh which calls QTLtools cis and creates permutations.txt files for different numbers of PC we need to compute how many significant cis QTLs are discovered for each number of PC

# column 19 i  permutation files: The P-value of association adjusted for the number of variants tested in cis given by the fitted beta distribution. We strongly recommend to use this adjusted P-value in any downstream analysis

CIS_QTL_FOLDER=/data/unige/funpopgen/davalos/project/GEUVADIS/step3_cisQTL
mkdir -p "$CIS_QTL_FOLDER/significant_cisQTLs" #create folder if not yet existing
SUMMARY_CIS_QTL_FILE=$CIS_QTL_FOLDER/significant_cisQTLs/summary.txt
echo -e "number_PC\tcis_QTL_found" > $SUMMARY_CIS_QTL_FILE # empty SUMMARY_CIS_QTL_FILE
PVALUE_THRESHOLD=0.05

for PC in 1 2 5 10 20 30 40 50 60 70 80 90 100; do

    FILE=$CIS_QTL_FOLDER/permutation_$PC.txt
    SIGNIFICANT_CISQTL_FILE=$CIS_QTL_FOLDER/significant_cisQTLs/sign_cisQTL_$PC.txt
    echo >$SIGNIFICANT_CISQTL_FILE # make sure the file is empty 
    COUNT_CIS_QTLS=0

    while read -r line; do 
        PVALUE=$(echo "$( echo "$line "| tr ' ' '\t' | cut -f19)")
        # echo "$PVALUE"

        if [[ "$PVALUE"<"$PVALUE_THRESHOLD" ]]; then
            COUNT_CIS_QTLS=$(expr $COUNT_CIS_QTLS + 1)
            echo $line>>$SIGNIFICANT_CISQTL_FILE    
        fi

    done < "$FILE"
    
    echo -e "$PC\t$COUNT_CIS_QTLS">>$SUMMARY_CIS_QTL_FILE

done


