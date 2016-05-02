#!/bin/bash

# Extract CHR/POS format from GoNL DNMs raw data
cut -f 2-3 /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs.txt | tail -n +2 | sort -k1,1n > /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs2.txt


# Combine per-chr DNMs (already formatted as CHR POS MU OBS)
cat /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs2_chr*_mu.txt > /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs2_mu.txt

# Combine per-chr subset data (already formatted as CHR POS MU OBS)
cat /net/bipolar/jedidiah/mutation/output/predicted/full/chr*_sub.txt > /net/bipolar/jedidiah/mutation/output/predicted/full/all_sub.txt

# grep -Fwf  /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs2.txt /net/bipolar/jedidiah/mutation/output/predicted/full/all_sub.txt | sed "s/$/\t1/" > /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs2_mu.txt

# Combine subset data with DNMs
cat /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs2_mu.txt /net/bipolar/jedidiah/mutation/output/predicted/full/all_sub.txt > /net/bipolar/jedidiah/mutation/output/predicted/full/rocdat_comb.txt
