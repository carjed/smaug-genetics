#!/bin/bash

chr=$1

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f1,3 chr${chr}_AT_CG_predicted.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_AT_CG.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f1,3 chr${chr}_AT_GC_predicted.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_AT_GC.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f1,3 chr${chr}_AT_TA_predicted.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_AT_TA.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f1,3 chr${chr}_GC_AT_predicted.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_GC_AT.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f1,3 chr${chr}_GC_CG_predicted.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_GC_CG.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f1,3 chr${chr}_GC_TA_predicted.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_GC_TA.wig

# echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f1,3 chr${chr}_GC_AT_predicted.txt) | /net/bipolar/jedidiah/mutation/output/predicted/wigToBigWig  https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes chr${chr}_GC_AT.bw
