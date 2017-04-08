#!/usr/local/bin/perl

##############################################################################
# Used to obtain full data from logistic regression model
# loops through reference genome and outputs 1 line per base, as long as
# valid covariate data exists
##############################################################################

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils 'pairwise';
use Cwd;
use Benchmark;
use Tie::File;
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

my $adj = 3;
my $data = $config->{data};
my $email = $config->{email};
my $parentdir = $config->{parentdir};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

my $today = POSIX::strftime('%Y%m%d', localtime);
my $slurmdir = "$parentdir/output/slurm/$today";
  make_path("$slurmdir");

my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );

foreach my $categ (@categs){
  my $fullfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_full.txt.gz";
  my $subcmd = "zcat $fullfile | awk '{print >> \"$parentdir/output/logmod_data/motifs/$categ/${categ}_\" \$4 \".txt\"}' &";

  forkExecWait($slurmcmd);
}


# forkExecWait($subcmd);
