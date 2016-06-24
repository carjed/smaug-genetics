#!/bin/bash

for i in /net/bipolar/bridges/sequencedata/umake/glfs/samples/chr1/1.5000000/*.glf;
do
  samtools-hybrid glfview $i  | cut -f1,2,4 > /net/bipolar/jedidiah/mutation/output/glf_depth/$i.dp
done
