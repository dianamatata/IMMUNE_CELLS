#!/bin/bash

for c in $(seq 1 22); do
	LINP=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/LCL_ALL.chr$c\.txt.gz
	FINP=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/FIB_ALL.chr$c\.txt.gz
	SINP=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/SUB_ALL.chr$c\.txt.gz

	LDIST=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/LCL_ALL.chr$c\.distance.txt.gz
	FDIST=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/FIB_ALL.chr$c\.distance.txt.gz
	SDIST=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/SUB_ALL.chr$c\.distance.txt.gz

	echo "zcat $LINP | awk '{ p1=(\$3+\$4)/2; p2=(\$8+\$9)/2; d=p2-p1; d=(d>0)?d:-d; if (d < 1000000) print \$0 }' | gzip -c > $LDIST" | bsub -q priority
	echo "zcat $FINP | awk '{ p1=(\$3+\$4)/2; p2=(\$8+\$9)/2; d=p2-p1; d=(d>0)?d:-d; if (d < 1000000) print \$0 }' | gzip -c > $FDIST" | bsub -q priority
	echo "zcat $SINP | awk '{ p1=(\$3+\$4)/2; p2=(\$8+\$9)/2; d=p2-p1; d=(d>0)?d:-d; if (d < 1000000) print \$0 }' | gzip -c > $SDIST" | bsub -q priority

	LSUB=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/LCL_ALL.chr$c\.subset.txt.gz
	FSUB=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/FIB_ALL.chr$c\.subset.txt.gz
	SSUB=/data/unige/funpopgen/odelanea/SGX/V2/data_correlation/cis/chromatin/SUB_ALL.chr$c\.subset.txt.gz

	echo "zcat $LINP | awk '{ d=\$6-\$1; if (d <= 250) print \$0 }' | gzip -c > $LSUB" | bsub -q priority
	echo "zcat $FINP | awk '{ d=\$6-\$1; if (d <= 250) print \$0 }' | gzip -c > $FSUB" | bsub -q priority
	echo "zcat $SINP | awk '{ d=\$6-\$1; if (d <= 250) print \$0 }' | gzip -c > $SSUB" | bsub -q priority
done

