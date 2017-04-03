#!/bin/bash

##############################################################################
# Script downloads and formats reference data
##############################################################################

#url=$1
#wget $url

refdir="/net/bipolar/jedidiah/mutation/reference_data"

# Reference genome
#wget -P $parentdir/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
wget -qO- ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz | gunzip > "$refdir/human_g1k_v37.fasta"

# hg19 chromosome lengths
wget -qO- https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes > "$refdir/hg19.genome"

# CpG islands
wget -qO- http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt | awk 'NR>1' | sort -k1,1 -k2,2n > "$refdir/cpg_islands_sorted.txt"

# Lamin-associated domains
wget -qO- http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/laminB1Lads.txt.gz | gunzip | awk 'NR>1 {print $2"\t"$3"\t"$4}' | bedtools sort -i - > "$refdir/lamin_B1_LADS2.bed"

# Histone marks
wget -r -nd -P . --accept-regex 'E062' http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/

gunzip *.broadPeak.gz
for f in *.broadPeak; do
	mv -- "$f" "${f%.broadPeak}.bed"
done

for i in E062*.bed; do
	bedtools sort -i $i > sort.$i
done

# DNase hypersensitive sites
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz | gunzip | cut -f1-3 | bedtools sort -i - > "$refdir/DHS.bed"

# Replication timing
wget -qO- http://mccarrolllab.com/wp-content/uploads/2015/03/Koren-et-al-Table-S2.zip | gunzip > "$refdir/lymph_rep_time.txt"

# Recombination rate
wget http://hgdownload.cse.ucsc.edu/gbdb/hg19/decode/SexAveraged.bw
bigWigToWig SexAveraged.bw SexAveraged.wig
echo "CHR\tSTART\tEND\tRATE" > recomb_rate.bed
awk 'NR>1' SexAveraged.wig | cat >> recomb_rate.bed

# GoNL de novo mutations
wget https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release5.2/GoNL_DNMs.txt > "$refdir/DNMs/GoNL_DNMs.txt"

# RefSeq exons