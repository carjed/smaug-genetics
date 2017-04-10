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

# hg19 chromosome lengths
curl -s "https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes" > "$refdir/hg19.genome"

# 1000G strict mask
curl -s  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed" | bedtools complement -i - -g "$refdir/hg19.genome" | bedtools sort | awk 'match($1, /chr[0-9]+$/) {print $0}' > "$refdir/testmask2.bed"

# Reference genome fasta
curl -s "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz" > "$refdir/human_g1k_v37/human_g1k_v37.fasta.gz"

for i in `seq 1 22`; do
  samtools faidx "$refdir/human_g1k_v37/human_g1k_v37.fasta.gz" $i | bgzip -c > "$refdir/human_g1k_v37/chr$i.fasta.gz"
  samtools faidx "$refdir/human_g1k_v37/chr$i.fasta.gz"
done

# Ancestral genome fasta
curl -s "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis_results/supporting/ancestral_alignments/human_ancestor_GRCh37_e59.tar.bz2" > "$refdir/human_ancestor_GRCh37_e59.tar.bz2"

tar -vjxf "$refdir/human_ancestor_GRCh37_e59.tar.bz2"

# ancdir="$refdir/human_ancestor_GRCh37_e59"
# mkdir $ancdir
# cd $ancdir
# tar -xvjf
# cd $refdir

# for i in `seq 1 22`; do
#   curl -s "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/technical/reference/ancestral_alignments/human_ancestor_$i.fa.bz2" | bzcat | sed "s,^>.*,>$CHR," | gzip -c > $OUT $ancdir/human_ancestor_$i.fa.gz
# done
#
# ls human_ancestor_*.fa.bz2 | while read IN; do
#     OUT=`echo $IN | sed 's,bz2$,gz,'`
#     CHR=`echo $IN | sed 's,human_ancestor_,, ; s,.fa.bz2,,'`
#     bzcat $IN | sed "s,^>.*,>$CHR," | gzip -c > $OUT
#     samtools faidx $OUT
# done

# mask fasta
bedtools maskfasta -fi "$refdir/human_g1k_v37.fasta" -bed "$refdir/testmask2.bed" -fo "$refdir/human_g1k_v37.premask.fasta"

perl -ane 'if(/\>/){$a++;print ">$a dna:chromosome\n"}else{print;}' "$refdir/human_g1k_v37.mask.fasta" > "$refdir/human_g1k_v37.mask.fasta"

rm -f "$refdir/human_g1k_v37.premask.fasta"

# CpG islands
curl -s  "http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt" | awk 'NR>1' | sort -k1,1 -k2,2n > "$refdir/cpg_islands_sorted.bed"

# Lamin-associated domains
curl -s  "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/laminB1Lads.txt.gz" | gunzip | awk 'NR>1 {print $2"\t"$3"\t"$4}' | bedtools sort -i - > "$refdir/lamin_B1_LADS2.bed"

# DNase hypersensitive sites
curl -s "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClusteredV3.bed.gz" | gunzip | cut -f1-3 | bedtools sort -i - > "$refdir/DHS.bed"

# Replication timing
curl -s "http://mccarrolllab.com/wp-content/uploads/2015/03/Koren-et-al-Table-S2.zip" | gunzip > "$refdir/lymph_rep_time.txt"

# Recombination rate
curl -s "http://hgdownload.cse.ucsc.edu/gbdb/hg19/decode/SexAveraged.bw" | "$refdir/SexAvaraged.bw"
bigWigToWig "$refdir/SexAveraged.bw" "$refdir/SexAveraged.wig"
echo "CHR\tSTART\tEND\tRATE" > "$refdir/recomb_rate.bed"
awk 'NR>1' "$refdir/SexAveraged.wig" | cat >> "$refdir/recomb_rate.bed"

# GoNL de novo mutations
curl -s "https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release5.2/GoNL_DNMs.txt" > "$refdir/DNMs/GoNL_DNMs.txt"

# ITMI de novo mutations
curl -s "http://www.nature.com/ng/journal/v48/n8/extref/ng.3597-S3.xlsx" > "$refdir/DNMs/goldmann_2016_dnms.xlsx"

# RefSeq v69 exons
# originally downloaded via UCSC genome browser;
# this command directly downloads the file used in analyses
curl -s "http://mutation.sph.umich.edu/hg19/GRCh37_RefSeq_sorted.bed" >  "$refdir/GRCh37_RefSeq_sorted.bed"

# Histone marks
wget -r -nd -P . --accept-regex 'E062' http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/

gunzip *.broadPeak.gz
for f in *.broadPeak; do
	mv -- "$f" "${f%.broadPeak}.bed"
done

for i in E062*.bed; do
	bedtools sort -i $i > sort.$i
done

# cytobands
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip > "$refdir/cytoBand.txt"

# get mask pct per band; used in downstream analysis
bedtools coverage -a "$refdir/testmask2.bed" -b "$refdir/cytoBand.txt" > "$refdir/testcov.bed"

# Human Accelerated Regions
curl -s "http://www.broadinstitute.org/ftp/pub/assemblies/mammals/29mammals/2xHARs.bed" > "$refdir/2xHARs.bed"

"$refdir/liftOver" "$refdir/2xHARs.bed" "$refdir/hg18ToHg19.over.chain.gz" "$refdir/2xHARs.hg19.bed" "$refdir/unlifted.bed"

bedtools sort -i "$refdir/2xHARs.hg19.bed" > "$refdir/2xHARs.hg19.sort.bed"

# Aggarwala & Voight rates
curl -s "http://www.nature.com/ng/journal/v48/n4/extref/ng.3511-S2.xlsx" > "$refdir/AV_rates.xlsx"

cd $curdir
