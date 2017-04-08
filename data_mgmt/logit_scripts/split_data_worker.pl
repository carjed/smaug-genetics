#!/usr/local/bin/perl

##############################################################################
# Split logit data by motif
##############################################################################
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

my $parentdir = $config->{parentdir};

my $catind = $ARGV[0]-1;
my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );
my $categ = $categs[$catind];

foreach my $chr (1..22){

  my $fullfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_full.txt.gz";
  my $subcmd = "zcat $fullfile | awk '{print >> \"$parentdir/output/logmod_data/motifs/$categ/${categ}_\" \$4 \".txt\"}'";

  forkExecWait($subcmd);
}