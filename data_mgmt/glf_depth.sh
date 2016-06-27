#!/bin/bash

##############################################################################
# Script scans glf files, and extracts only the relevant info for every
# 10th base
##############################################################################

for i in /net/bipolar/bridges/sequencedata/umake/glfs/samples/chr1/1.5000000/*.glf;
do
  fname=$(basename $i)
  samtools-hybrid glfview $i  | cut -f1-4 | awk '$2%10==0 && $3 ~ /[ACGT]/' > /net/bipolar/jedidiah/mutation/output/glf_depth/$fname.dp
done
