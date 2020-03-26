#!/bin/bash

FOLDER=/home/users/a/avalosma/scratch/Gene_CRD/quantify

for cell_type_data in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
        for cell_type_crd in 'EGAD00001002670'  'EGAD00001002672' 'EGAD00001002673' ; do
                if [[ "$cell_type_data" != "$cell_type_crd" ]] ; then
                        name1=$(echo ${cell_type_data: (-2)})data_vs_$(echo ${cell_type_crd: (-2)})crd
                        name=$FOLDER/$name1
			echo $name
			gunzip ${name}.chr*.mean.txt.gz
			awk 'FNR>1||NR==1' ${name}.chr*.mean.txt > ${name}.ALLchr.mean.unsorted.txt
			awk 'NR == 1; NR > 1 {print $0 | "sort -V -k1,1 -k2,2n"}' ${name}.ALLchr.mean.unsorted.txt | bgzip -c >	${name}.ALLchr.mean.txt.gz	
			tabix -p bed ${name}.ALLchr.mean.txt.gz
			bgzip ${name}.chr*.mean.txt
			rm ${name}.ALLchr.mean.unsorted.txt
		fi
	done
done

