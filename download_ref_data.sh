#!/bin/bash

#############################################################################
# Script downloads and formats reference data
#############################################################################

# Function from https://gist.github.com/pkuczynski/8665367
parse_yaml() {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\)\($w\)$s:$s\"\(.*\)\"$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

# read yaml file
eval $(parse_yaml _config.yaml "config_")

refdir="$config_parentdir/reference_data"
mkdir $refdir

curdir=${PWD}
# echo $curdir

cd $refdir

#############################################################################
# hg19 chromosome lengths
#############################################################################
curl -s "https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes" > "$refdir/hg19.genome"

# Make fixed-width windows
bedtools makewindows -g "$refdir/hg19.genome" -w 1000000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "$refdir/genome.1000kb.sorted.bed"

bedtools makewindows -g "$refdir/hg19.genome" -w 5000000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "$refdir/genome.5000kb.sorted.bed"

bedtools makewindows -g "$refdir/hg19.genome" -w 100000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "$refdir/genome.100kb.sorted.bed"

bedtools makewindows -g "$refdir/hg19.genome" -w 10000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "$refdir/genome.10kb.sorted.bed"

bedtools makewindows -g "$refdir/hg19.genome" -w 3000000000 | grep -Ev "_|X|Y|M" | sort -k 1,1 -k2,2n > "$refdir/genome.full.sorted.bed"

#############################################################################
# 1000G strict mask
#############################################################################
curl -s  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed" | bedtools complement -i - -g "$refdir/hg19.genome" | bedtools sort | awk 'match($1, /chr[0-9]+$/) {print $0}' > "$refdir/testmask2.bed"

#############################################################################
# Reference genomes in fasta format
#############################################################################
# v37
mkdir "$refdir/human_g1k_v37"
curl -s "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz" | gunzip -c > "$refdir/human_g1k_v37/human_g1k_v37.fasta"

for i in `seq 1 22`; do
  samtools faidx "$refdir/human_g1k_v37/human_g1k_v37.fasta" $i | bgzip -c > "$refdir/human_g1k_v37/chr$i.fasta.gz"
  samtools faidx "$refdir/human_g1k_v37/chr$i.fasta.gz"
done

# mask v37
mkdir "$refdir/human_g1k_v37_mask"
bedtools maskfasta -fi "$refdir/human_g1k_v37/human_g1k_v37.fasta" -bed "$refdir/testmask2.bed" -fo "$refdir/human_g1k_v37_mask/human_g1k_v37.mask.fasta"

# perl -ane 'if(/\>/){$a++;print ">$a dna:chromosome\n"}else{print;}' "$refdir/human_g1k_v37_mask/human_g1k_v37.mask.fasta" > "$refdir/human_g1k_v37_mask/human_g1k_v37.mask.fasta"
#
# rm -f "$refdir/human_g1k_v37_mask/human_g1k_v37.premask.fasta"

for i in `seq 1 22`; do
  samtools faidx "$refdir/human_g1k_v37_mask/human_g1k_v37.mask.fasta" $i | bgzip -c > "$refdir/human_g1k_v37_mask/chr$i.fasta.gz"
  samtools faidx "$refdir/human_g1k_v37_mask/chr$i.fasta.gz"
done

# Ancestral genome
curl -s "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2" > "$refdir/human_ancestor_GRCh37_e59.tar.bz2"

tar -vjxf "$refdir/human_ancestor_GRCh37_e59.tar.bz2"

#############################################################################
# GC content in 10kb windows
#############################################################################
bedtools nuc -fi "$refdir/human_g1k_v37/human_g1k_v37.fasta" -bed <(sed s/chr// "$refdir/genome.10kb.sorted.bed") > "$refdir/gc10kb.bed"

#############################################################################
# CpG islands
#############################################################################
curl -s  "http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt" | awk 'NR>1' | sort -k1,1 -k2,2n > "$refdir/cpg_islands_sorted.bed"

#############################################################################
# Lamin-associated domains
#############################################################################
curl -s  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/laminB1Lads.txt.gz" | gunzip | awk 'NR>1 {print $2"\t"$3"\t"$4}' | bedtools sort -i - > "$refdir/lamin_B1_LADS2.bed"

#############################################################################
# DNase hypersensitive sites
#############################################################################
curl -s "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz" | gunzip | cut -f1-3 | bedtools sort -i - > "$refdir/DHS.bed"

#############################################################################
# Replication timing
#############################################################################
curl -s "http://mccarrolllab.com/wp-content/uploads/2015/03/Koren-et-al-Table-S2.zip" | gunzip > "$refdir/lymph_rep_time.txt"

#############################################################################
# Recombination rate
#############################################################################
curl -s "http://hgdownload.cse.ucsc.edu/gbdb/hg19/decode/SexAveraged.bw" > "$refdir/SexAveraged.bw"
bigWigToWig "$refdir/SexAveraged.bw" "$refdir/SexAveraged.wig"
echo "CHR\tSTART\tEND\tRATE" > "$refdir/recomb_rate.bed"
awk 'NR>1' "$refdir/SexAveraged.wig" | cat >> "$refdir/recomb_rate.bed"

#############################################################################
# Histone marks
#############################################################################
wget -r -nd -P . --accept-regex 'E062' https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/

gunzip *.broadPeak.gz
for f in *.broadPeak; do
	mv -- "$f" "${f%.broadPeak}.bed"
done

for i in E062*.bed; do
	bedtools sort -i $i > sort.$i
done

#############################################################################
# de novo mutations
#############################################################################
mkdir $refdir/DNMs
# GoNL
curl -s "https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release5.2/GoNL_DNMs.txt" > "$refdir/DNMs/GoNL_DNMs.txt"

# ITMI
curl -s "http://www.nature.com/ng/journal/v48/n8/extref/ng.3597-S3.xlsx" > "$refdir/DNMs/goldmann_2016_dnms.xlsx"

#############################################################################
# RefSeq v69 exons
# originally downloaded via UCSC genome browser;
# this command directly downloads the file used in analyses
#############################################################################
curl -s "http://mutation.sph.umich.edu/hg19/GRCh37_RefSeq_sorted.bed" >  "$refdir/GRCh37_RefSeq_sorted.bed"

#############################################################################
# cytobands
#############################################################################
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip > "$refdir/cytoBand.txt"

# get mask pct per band; used in downstream analysis
bedtools coverage -a "$refdir/testmask2.bed" -b "$refdir/cytoBand.txt" > "$refdir/testcov.bed"

#############################################################################
# Human Accelerated Regions
#############################################################################
curl -s "ftp://ftp.broadinstitute.org/pub/assemblies/mammals/29mammals/2xHARs.bed" > "$refdir/2xHARs.bed"

curl -s "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver" > liftOver
chmod +x liftOver

curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz" > "$refdir/hg18ToHg19.over.chain.gz"

"$refdir/liftOver" "$refdir/2xHARs.bed" "$refdir/hg18ToHg19.over.chain.gz" "$refdir/2xHARs.hg19.bed" "$refdir/unlifted.bed"

bedtools sort -i "$refdir/2xHARs.hg19.bed" > "$refdir/2xHARs.hg19.sort.bed"

#############################################################################
# Aggarwala & Voight rates
#############################################################################
curl -s "http://www.nature.com/ng/journal/v48/n4/extref/ng.3511-S2.xlsx" > "$refdir/AV_rates.xlsx"

cd $curdir
