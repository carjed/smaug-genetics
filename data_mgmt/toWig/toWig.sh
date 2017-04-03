#!/bin/bash

chr=$1

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}.AT_CG.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_AT_CG.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}.AT_GC.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_AT_GC.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}.AT_TA.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_AT_TA.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}.GC_AT.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_GC_AT.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}.GC_CG.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_GC_CG.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}.GC_TA.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_GC_TA.wig

# echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 chr${chr}_GC_AT_predicted.txt) | /net/bipolar/jedidiah/mutation/output/predicted/wigToBigWig  https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes chr${chr}_GC_AT.bw
