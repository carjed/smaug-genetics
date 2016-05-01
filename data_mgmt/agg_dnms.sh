#!/bin/bash

cat /net/bipolar/jedidiah/mutation/output/predicted/full/chr*_sub.txt > /net/bipolar/jedidiah/mutation/output/predicted/full/all_sub.txt

grep -Fwf  /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs2.txt /net/bipolar/jedidiah/mutation/output/predicted/full/all_sub.txt | sed "s/$/\t1/" > GoNL_DNMs2_mu.txt

cat /net/bipolar/jedidiah/mutation/reference_data/GoNL_DNMs2_mu.txt /net/bipolar/jedidiah/mutation/output/predicted/full/all_sub.txt > /net/bipolar/jedidiah/mutation/output/predicted/full/rocdat_comb.txt
