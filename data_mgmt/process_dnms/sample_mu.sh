#!/bin/bash

##############################################################################
# Step 0:
# Merges per-category predicted rates into single chromosome file, sorted by
# position
##############################################################################

categ=$1

# Get random selection of sites with predicted rates
# Must modify to include category
awk -v categ="$categ" 'BEGIN {srand()} !/^$/ { if (rand() <= .005) print chr"\t"$0"\t"0"\t"categ}' /net/bipolar/jedidiah/mutation/output/predicted/chr*.${categ}.txt > /net/bipolar/jedidiah/mutation/output/predicted/${categ}.sub.txt
