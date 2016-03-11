#!/bin/bash

chr=$1

paste /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_AT_CG_predicted.txt <(cut -d $'\t' -f3 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_AT_GC_predicted.txt) > /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_tmp2.txt
paste /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_tmp2.txt <(cut -d $'\t' -f3 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_AT_TA_predicted.txt) | awk -v OFS='\t' '{print $1,$3+$4+$5}' > /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_AT.txt

paste /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_GC_AT_predicted.txt <(cut -d $'\t' -f3 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_GC_CG_predicted.txt) > /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_tmp.txt
paste /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_tmp.txt <(cut -d $'\t' -f3 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_GC_TA_predicted.txt) | awk -v OFS='\t' '{print $1,$3+$4+$5}' > /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_GC.txt

cat /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_AT.txt /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_GC.txt | sort -n -k1,1 > /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_full.txt

rm -f /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_{AT,GC,tmp,tmp2}.txt

awk -F$'\t' '{print ${chr}"\t"($1-1)"\t"$1"\t"$2}' /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_full.txt > /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_full.bed

bedtools intersect -a /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_full.bed -b /net/bipolar/jedidiah/mutation/reference_data/20140520.strict_mask2.bed  | awk '{print "\t"$3"\t"$4}' > /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_full_mask.txt
