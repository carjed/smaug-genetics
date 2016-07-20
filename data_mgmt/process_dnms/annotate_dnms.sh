#!/bin/bash

##############################################################################
# Step 0:
# Merges per-category predicted rates into single chromosome file, sorted by
# position
##############################################################################

# chr=$1
categ=$1

for i in {1:22}
do
  grep -Fwf  <(grep "\s${i}\s" /net/bipolar/jedidiah/mutation/reference_data/DNMs/GoNL_${categ}.txt | cut -f 3)  /net/bipolar/jedidiah/mutation/output/predicted/chr${i}.${categ}.txt | awk -v categ="$categ" '{print $0"\t"1"\tcateg"}' >> /net/bipolar/jedidiah/mutation/reference_data/DNMs/GoNL_${categ}.anno.txt
done
