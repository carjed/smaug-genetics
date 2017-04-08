#!/usr/local/bin/perl

use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Path qw(make_path);
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));

my $config = LoadFile("$configpath/_config.yaml");

my $adj = $config->{adj};
my $mac = $config->{mac};
my $binw = $config->{binw};
my $data = $config->{data};
my $bin_scheme = $config->{bin_scheme};
my $parentdir = $config->{parentdir};
my $count_motifs = $config->{count_motifs};
my $expand_summ = $config->{expand_summ};

my $chr=$ARGV[0];

my $atcgcmd = "echo -e \"variableStep\tchrom=chr${chr}\" | cat - <(cut -f2,4 $parentdir/output/predicted/chr${chr}.AT_CG.txt) > $parentdir/output/predicted/tracks/chr${chr}_AT_CG.wig";

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}.AT_GC.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_AT_GC.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}.AT_TA.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_AT_TA.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}.GC_AT.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_GC_AT.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}.GC_CG.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_GC_CG.wig

echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 /net/bipolar/jedidiah/mutation/output/predicted/chr${chr}.GC_TA.txt) > /net/bipolar/jedidiah/mutation/output/predicted/tracks/chr${chr}_GC_TA.wig

# echo -e "variableStep\tchrom=chr${chr}" | cat - <(cut -f2,4 chr${chr}_GC_AT_predicted.txt) | /net/bipolar/jedidiah/mutation/output/predicted/wigToBigWig  https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes chr${chr}_GC_AT.bw
