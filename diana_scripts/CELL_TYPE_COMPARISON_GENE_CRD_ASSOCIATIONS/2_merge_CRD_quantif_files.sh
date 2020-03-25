# gzip -d 70_vs_72_ALL.chr*.mean.txt.gz
# awk 'FNR>1||NR==1' 70_vs_72_ALL.chr*.mean.txt > 70_vs_72_ALL.ALLchr.mean.unsorted.txt
# (head -1 70_vs_72_ALL.ALLchr.mean.unsorted.txt && tail -n +2 70_vs_72_ALL.ALLchr.mean.unsorted.txt | sort -V -k1,1 -k2,2n) | bgzip -c > 70_vs_72_ALL.ALLchr.mean.txt.gz
# tabix -p bed 70_vs_72_ALL.ALLchr.mean.txt.gz
# gzip 70_vs_72_ALL.chr*.mean.txt
#
# gzip -d 70_vs_73_ALL.chr*.mean.txt.gz
# awk 'FNR>1||NR==1' 70_vs_73_ALL.chr*.mean.txt > 70_vs_73_ALL.ALLchr.mean.unsorted.txt
# (head -1 70_vs_73_ALL.ALLchr.mean.unsorted.txt && tail -n +2 70_vs_73_ALL.ALLchr.mean.unsorted.txt | sort -V -k1,1 -k2,2n) | bgzip -c > 70_vs_73_ALL.ALLchr.mean.txt.gz
# tabix -p bed 70_vs_73_ALL.ALLchr.mean.txt.gz
# gzip 70_vs_73_ALL.chr*.mean.txt

gzip -d 72_vs_73_ALL.chr*.mean.txt.gz
awk 'FNR>1||NR==1' 72_vs_73_ALL.chr*.mean.txt > 72_vs_73_ALL.ALLchr.mean.unsorted.txt
(head -1 72_vs_73_ALL.ALLchr.mean.unsorted.txt && tail -n +2 72_vs_73_ALL.ALLchr.mean.unsorted.txt | sort -V -k1,1 -k2,2n) | bgzip -c > 72_vs_73_ALL.ALLchr.mean.txt.gz
tabix -p bed 72_vs_73_ALL.ALLchr.mean.txt.gz
gzip 72_vs_73_ALL.chr*.mean.txt

gzip -d 73_vs_72_ALL.chr*.mean.txt.gz
awk 'FNR>1||NR==1' 73_vs_72_ALL.chr*.mean.txt > 73_vs_72_ALL.ALLchr.mean.unsorted.txt
(head -1 73_vs_72_ALL.ALLchr.mean.unsorted.txt && tail -n +2 73_vs_72_ALL.ALLchr.mean.unsorted.txt | sort -V -k1,1 -k2,2n) | bgzip -c > 73_vs_72_ALL.ALLchr.mean.txt.gz
tabix -p bed 73_vs_72_ALL.ALLchr.mean.txt.gz
gzip 73_vs_72_ALL.chr*.mean.txt

gzip -d 72_vs_70_ALL.chr*.mean.txt.gz
awk 'FNR>1||NR==1' 72_vs_70_ALL.chr*.mean.txt > 72_vs_70_ALL.ALLchr.mean.unsorted.txt
(head -1 72_vs_70_ALL.ALLchr.mean.unsorted.txt && tail -n +2 72_vs_70_ALL.ALLchr.mean.unsorted.txt | sort -V -k1,1 -k2,2n) | bgzip -c > 72_vs_70_ALL.ALLchr.mean.txt.gz
tabix -p bed 72_vs_70_ALL.ALLchr.mean.txt.gz
gzip 72_vs_70_ALL.chr*.mean.txt

gzip -d 73_vs_70_ALL.chr*.mean.txt.gz
awk 'FNR>1||NR==1' 73_vs_70_ALL.chr*.mean.txt > 73_vs_70_ALL.ALLchr.mean.unsorted.txt
(head -1 73_vs_70_ALL.ALLchr.mean.unsorted.txt && tail -n +2 73_vs_70_ALL.ALLchr.mean.unsorted.txt | sort -V -k1,1 -k2,2n) | bgzip -c > 73_vs_70_ALL.ALLchr.mean.txt.gz
tabix -p bed 73_vs_70_ALL.ALLchr.mean.txt.gz
gzip 73_vs_70_ALL.chr*.mean.txt
