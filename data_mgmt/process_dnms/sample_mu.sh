#!/bin/bash

##############################################################################
# Step 0:
# Merges per-category predicted rates into single chromosome file, sorted by
# position
##############################################################################

categ=$1

# Get random selection of sites with predicted rates
# Must modify to include category
awk -v categ="$categ" 'BEGIN {srand()} !/^$/ { if (rand() <= .005) print $0"\t"0"\t"categ}' /net/bipolar/jedidiah/mutation/output/predicted/chr*.${categ}.txt > /net/bipolar/jedidiah/mutation/output/predicted/${categ}.sub.txt

# Combine null sites with DNMs
cat /net/bipolar/jedidiah/mutation/output/predicted/*.sub.txt /net/bipolar/jedidiah/mutation/reference_data/DNMs/*.anno.txt | sort -V -k1,2 > /net/bipolar/jedidiah/mutation/output/rocdat.sort.txt

# Annotate with rate tables (can omit and do this in R after selecting 1M sites)
# for faster processing
# Also need to specify how {group}_7bp_rates.txt tables are generated
# CMD to annotate with common rates
# perl data_mgmt/process_dnms/anno_rate.pl --adj 3 --in /net/bipolar/jedidiah/mutation/output/rocdat.sort.txt --rates /net/bipolar/jedidiah/mutation/common_7bp_rates.txt --seq --out /net/bipolar/jedidiah/mutation/output/rocdat.7bp.1.txt

# Annotate with ERV rates
# perl data_mgmt/process_dnms/anno_rate.pl --adj 3 --in /net/bipolar/jedidiah/mutation/output/rocdat.7bp.1.txt --rates /net/bipolar/jedidiah/mutation/ERV_7bp_rates.txt --out /net/bipolar/jedidiah/mutation/output/rocdat.7bp.txt

# Annotate with AV rates
# perl data_mgmt/process_dnms/anno_rate.pl --adj 3 --in /net/bipolar/jedidiah/mutation/output/rocdat.7bp.1.txt --rates /net/bipolar/jedidiah/mutation/av_7bp_rates.txt --out /net/bipolar/jedidiah/mutation/output/rocdat.7bp.2.txt
