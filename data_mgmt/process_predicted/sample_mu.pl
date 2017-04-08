#!/usr/local/bin/perl

##############################################################################
# Step 0:
# Sample sites for non-mutated background
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

my $adj = $config->{adj};
my $mac = $config->{mac};
my $binw = $config->{binw};
my $data = $config->{data};
my $bin_scheme = $config->{bin_scheme};
my $parentdir = $config->{parentdir};
my $count_motifs = $config->{count_motifs};
my $expand_summ = $config->{expand_summ};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );

print "Preparing de novo data...\n";
my $prepdnmcmd = "Rscript $parentdir/smaug-genetics/R/read_dnms.r TRUE $parentdir";
forkExecWait($prepdnmcmd);

foreach my $categ (@categs) {
  # Get random selection of sites with predicted rates
  # Must modify to include category
  print "Subsetting $categ sites...\n";
  my $samplecmd = "awk -v categ=\"$categ\" 'BEGIN {srand()} !/^\$/ { if (rand() <= .005) print \$0\"\\t\"0\"\\t\"categ}' $parentdir/output/predicted/chr*.${categ}.txt > $parentdir/output/predicted/${categ}.sub_new.txt";
  forkExecWait($samplecmd);

  # remove per-category DNM data if it already exists
  my $cleanupcmd = "rm -f $parentdir/reference_data/DNMs/GoNL_${categ}.anno.txt";
  forkExecWait($cleanupcmd);

  # Testing: get de novo data in same script
  for my $i (1 .. 22) {
    my $tmpfile = "$parentdir/reference_data/DNMs/tmp.txt";

    my $buildquerycmd = "grep \"\\s$i\\s\" $parentdir/reference_data/DNMs/GoNL_${categ}.txt | cut -f 3 > $tmpfile";
    forkExecWait($buildquerycmd);

    my $dnmannocmd = "grep -Fwf  $tmpfile $parentdir/output/predicted/chr$i.${categ}.txt | awk -v categ=\"$categ\" '{print \$0\"\\t\"1\"\\t\"categ}' >> $parentdir/reference_data/DNMs/GoNL_${categ}.anno.txt";
    forkExecWait($dnmannocmd);
  }
}

print "Combining and sorting data...\n";
# Combine null sites with DNMs
my $sortcmd = "cat $parentdir/output/predicted/*.sub_new.txt $parentdir/reference_data/DNMs/*.anno.txt | sort -V -k1,2 > $parentdir/output/rocdat.sort_new.txt";
forkExecWait($sortcmd);
