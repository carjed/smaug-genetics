#!/bin/bash

paste /net/bipolar/jedidiah/mutation/output/predicted/chr$1_AT_CG_predicted.txt <(cut -d $'\t' -f3 /net/bipolar/jedidiah/mutation/output/predicted/chr$1_AT_GC_predicted.txt) > /net/bipolar/jedidiah/mutation/output/predicted/chr$1_tmp2.txt
paste /net/bipolar/jedidiah/mutation/output/predicted/chr$1_tmp2.txt <(cut -d $'\t' -f3 /net/bipolar/jedidiah/mutation/output/predicted/chr$1_AT_TA_predicted.txt) | awk -v OFS='\t' '{print $1,$3+$4+$5}' > /net/bipolar/jedidiah/mutation/output/predicted/chr$1_AT.txt

paste /net/bipolar/jedidiah/mutation/output/predicted/chr$1_GC_AT_predicted.txt <(cut -d $'\t' -f3 /net/bipolar/jedidiah/mutation/output/predicted/chr$1_GC_CG_predicted.txt) > /net/bipolar/jedidiah/mutation/output/predicted/chr$1_tmp.txt
paste /net/bipolar/jedidiah/mutation/output/predicted/chr$1_tmp.txt <(cut -d $'\t' -f3 /net/bipolar/jedidiah/mutation/output/predicted/chr$1_GC_TA_predicted.txt) | awk -v OFS='\t' '{print $1,$3+$4+$5}' > /net/bipolar/jedidiah/mutation/output/predicted/chr$1_GC.txt

cat /net/bipolar/jedidiah/mutation/output/predicted/chr$1_AT.txt /net/bipolar/jedidiah/mutation/output/predicted/chr$1_GC.txt | sort -n -k1,1 > /net/bipolar/jedidiah/mutation/output/predicted/chr$1_full.txt

rm -f /net/bipolar/jedidiah/mutation/output/predicted/chr$1_{AT,GC,tmp,tmp2}.txt

awk -F$'\t' '{print ($1-1)"\t"$1"\t"$2}' /net/bipolar/jedidiah/mutation/output/predicted/chr$1_full.txt > /net/bipolar/jedidiah/mutation/output/predicted/chr$1_full.bed

bedtools intersect -a /net/bipolar/jedidiah/mutation/output/predicted/chr$1_full.bed -b /net/bipolar/jedidiah/mutation/reference_data/20140520.strict_mask2.bed > /net/bipolar/jedidiah/mutation/output/predicted/chr$1_full_mask.txt
