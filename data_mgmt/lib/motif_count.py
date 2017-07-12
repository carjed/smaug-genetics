#!/usr/bin/python

from __future__ import print_function
import os
import sys
import textwrap
import argparse
import itertools
import timeit
import time
import csv
# import numpy as np
from subprocess import call
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
# import nimfa
# from util import *

###############################################################################
# Parse arguments
###############################################################################
# start = timeit.default_timer()

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input",
                    help="input fasta file",
                    required=True,
                    type=str,
                    default=sys.stdin)

parser.add_argument("-m", "--motifs",
                    help="input motif list",
                    required=True,
                    type=str)

parser.add_argument("-o", "--output",
                    help="output file",
                    required=True,
                    type=str)

args = parser.parse_args()

def occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

motif_dict = {}
with open(args.motifs) as f:
    motif_list = f.read().splitlines()

for m in motif_list:
    motif_dict[m] = 0

fasta_reader = Fasta(args.input, read_ahead=10000)

count = 0
for key in fasta_reader.keys():
    seq = fasta_reader[key]
    seqstr = seq[0:len(seq)].seq

    for m in motif_dict.keys():
        occ = occurrences(seqstr, m)
        motif_dict[m] += occ
    count += 1
    # if count % 1000 == 0:
    print("processed", count, "records in fasta")


outfile = open(args.output, 'w')
writer = csv.writer(outfile, delimiter = '\t')
for key, value in motif_dict.iteritems():
    writer.writerow([key] + [value])

# for key in motif_dict.keys():
#     outstr = key +
#     outfile.write()
# print(motif_dict)
