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

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait);

my $categ = $ARGV[0];
my $chr = $ARGV[1];

my $outpath = "$parentdir/output/logmod_data/motifs2/$categ";

my $fullfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_full.txt.gz";
my $subcmd = "zcat $fullfile | awk '{print >> \"$outpath/${categ}_\" substr(\$3, 1, 7) \".txt\"}'";
forkExecWait($subcmd);
