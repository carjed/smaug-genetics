#!/bin/bash

chr=$1

# cut -f 2-3 /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs.txt | tail -n +2 | sort -k1,1n  | sed "s/$/\t1/" | cut -f 1-3 > /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs2.txt

# | sed "s/$/\t0/"
# awk -v chr="$chr" 'BEGIN {srand()} !/^$/ { if (rand() <= .005) print chr"\t"$0"\t"0}' /net/bipolar/jedidiah/mutation/output/predicted/full/chr${chr}_full.txt > /net/bipolar/jedidiah/mutation/output/predicted/full/chr${chr}_sub.txt

grep -Fwf  <(grep "^$1" /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs2.txt | cut -f 2)  /net/bipolar/jedidiah/mutation/output/predicted/full/chr${chr}_full.txt | sed "s/$1\t$/\t1/" > /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs2_chr${chr}_mu.txt

# cat /net/bipolar/jedidiah/mutation/output/predicted/full/chr*_sub.txt > /net/bipolar/jedidiah/mutation/output/predicted/full/all_sub.txt
