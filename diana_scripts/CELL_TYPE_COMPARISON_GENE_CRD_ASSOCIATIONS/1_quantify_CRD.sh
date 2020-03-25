#!/bin/bash

BIN=/data/unige/funpopgen/odelanea/SHARE/clomics/bin/clomics

DATAFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0
CRDFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002672_CLOMICS_v3.0
OUTFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/

LI=$DATAFOLDER/merged_residuals.bed.gz

LM=$CRDFOLDER/EGAD00001002672_ALL.modules.MOD1.NRE2.txt.gz

for c in $(seq 1 22); do
        LT=$CRDFOLDER/EGAD00001002672_ALL.chr$c\.module.txt.gz

        LOP=$OUTFOLDER/70_vs_72_ALL.chr$c\.pc1
        LOM=$OUTFOLDER/70_vs_72_ALL.chr$c\.mean
        LOS=$OUTFOLDER/70_vs_72_ALL.chr$c\.loom

        #echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOP\.txt.gz --pca 1 --normal" | bsub -o $LOP\.log -R "rusage[mem=8000]" -M 8000000
        #echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOM\.txt.gz --mean --normal" | bsub -o $LOM\.log -R "rusage[mem=8000]" -M 8000000
        #echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOS\.txt.gz --loo 0 --normal" | bsub -o $LOS\.log -R "rusage[mem=8000]" -M 8000000

done

DATAFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0
CRDFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002673_CLOMICS_v3.0
OUTFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/

LI=$DATAFOLDER/merged_residuals.bed.gz

LM=$CRDFOLDER/EGAD00001002673_ALL.modules.MOD1.NRE2.txt.gz

for c in $(seq 1 22); do
        LT=$CRDFOLDER/EGAD00001002673_ALL.chr$c\.module.txt.gz

        LOP=$OUTFOLDER/70_vs_73_ALL.chr$c\.pc1
        LOM=$OUTFOLDER/70_vs_73_ALL.chr$c\.mean
        LOS=$OUTFOLDER/70_vs_73_ALL.chr$c\.loom

        #echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOP\.txt.gz --pca 1 --normal" | bsub -o $LOP\.log -R "rusage[mem=8000]" -M 8000000
        #echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOM\.txt.gz --mean --normal" | bsub -o $LOM\.log -R "rusage[mem=8000]" -M 8000000
        #echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOS\.txt.gz --loo 0 --normal" | bsub -o $LOS\.log -R "rusage[mem=8000]" -M 8000000

done

DATAFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002672_CLOMICS_v3.0
CRDFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002673_CLOMICS_v3.0
OUTFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/

LI=$DATAFOLDER/merged_residuals.bed.gz

LM=$CRDFOLDER/EGAD00001002673_ALL.modules.MOD1.NRE2.txt.gz

for c in $(seq 1 22); do
        LT=$CRDFOLDER/EGAD00001002673_ALL.chr$c\.module.txt.gz

        LOP=$OUTFOLDER/72_vs_73_ALL.chr$c\.pc1
        LOM=$OUTFOLDER/72_vs_73_ALL.chr$c\.mean
        LOS=$OUTFOLDER/72_vs_73_ALL.chr$c\.loom

        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOP\.txt.gz --pca 1 --normal" | bsub -o $LOP\.log -R "rusage[mem=8000]" -M 8000000
        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOM\.txt.gz --mean --normal" | bsub -o $LOM\.log -R "rusage[mem=8000]" -M 8000000
        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOS\.txt.gz --loo 0 --normal" | bsub -o $LOS\.log -R "rusage[mem=8000]" -M 8000000

done



DATAFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002673_CLOMICS_v3.0
CRDFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002672_CLOMICS_v3.0
OUTFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/

LI=$DATAFOLDER/merged_residuals.bed.gz

LM=$CRDFOLDER/EGAD00001002672_ALL.modules.MOD1.NRE2.txt.gz

for c in $(seq 1 22); do
        LT=$CRDFOLDER/EGAD00001002672_ALL.chr$c\.module.txt.gz

        LOP=$OUTFOLDER/73_vs_72_ALL.chr$c\.pc1
        LOM=$OUTFOLDER/73_vs_72_ALL.chr$c\.mean
        LOS=$OUTFOLDER/73_vs_72_ALL.chr$c\.loom

        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOP\.txt.gz --pca 1 --normal" | bsub -o $LOP\.log -R "rusage[mem=8000]" -M 8000000
        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOM\.txt.gz --mean --normal" | bsub -o $LOM\.log -R "rusage[mem=8000]" -M 8000000
        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOS\.txt.gz --loo 0 --normal" | bsub -o $LOS\.log -R "rusage[mem=8000]" -M 8000000

done

DATAFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002673_CLOMICS_v3.0
CRDFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0
OUTFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/

LI=$DATAFOLDER/merged_residuals.bed.gz

LM=$CRDFOLDER/EGAD00001002670_ALL.modules.MOD1.NRE2.txt.gz

for c in $(seq 1 22); do
        LT=$CRDFOLDER/EGAD00001002670_ALL.chr$c\.module.txt.gz

        LOP=$OUTFOLDER/73_vs_70_ALL.chr$c\.pc1
        LOM=$OUTFOLDER/73_vs_70_ALL.chr$c\.mean
        LOS=$OUTFOLDER/73_vs_70_ALL.chr$c\.loom

        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOP\.txt.gz --pca 1 --normal" | bsub -o $LOP\.log -R "rusage[mem=8000]" -M 8000000
        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOM\.txt.gz --mean --normal" | bsub -o $LOM\.log -R "rusage[mem=8000]" -M 8000000
        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOS\.txt.gz --loo 0 --normal" | bsub -o $LOS\.log -R "rusage[mem=8000]" -M 8000000

done


DATAFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002672_CLOMICS_v3.0
CRDFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/EGAD00001002670_CLOMICS_v3.0
OUTFOLDER=/data/unige/funpopgen/grey2/SYSCID/BLUEPRINT_DATA/CRD/THREE_CELL_TYPES/CLOMICS/CELL_TYPE_COMPARISON_GENE_CRD_ASSOCIATIONS/quantify/

LI=$DATAFOLDER/merged_residuals.bed.gz

LM=$CRDFOLDER/EGAD00001002670_ALL.modules.MOD1.NRE2.txt.gz

for c in $(seq 1 22); do
        LT=$CRDFOLDER/EGAD00001002670_ALL.chr$c\.module.txt.gz

        LOP=$OUTFOLDER/72_vs_70_ALL.chr$c\.pc1
        LOM=$OUTFOLDER/72_vs_70_ALL.chr$c\.mean
        LOS=$OUTFOLDER/72_vs_70_ALL.chr$c\.loom

        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOP\.txt.gz --pca 1 --normal" | bsub -o $LOP\.log -R "rusage[mem=8000]" -M 8000000
        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOM\.txt.gz --mean --normal" | bsub -o $LOM\.log -R "rusage[mem=8000]" -M 8000000
        echo "$BIN quantify --bed $LI --region $c --tree $LT $LM --out $LOS\.txt.gz --loo 0 --normal" | bsub -o $LOS\.log -R "rusage[mem=8000]" -M 8000000

done


