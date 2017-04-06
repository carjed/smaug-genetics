#!/usr/local/bin/perl

##############################################################################
# Step 0:
# Merges per-category predicted rates into single chromosome file, sorted by
# position
##############################################################################

##############################################################################
#Initialize inputs, options, and defaults and define errors
##############################################################################
use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Cwd;
use Benchmark;
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname($relpath);

my $config = LoadFile("$configpath/_config.yaml");

my $adj = $config->{adj};
my $mac = $config->{mac};
my $binw = $config->{binw};
my $data = $config->{data};
my $bin_scheme = $config->{bin_scheme};
my $parentdir = $config->{parentdir};
my $count_motifs = $config->{count_motifs};
my $expand_summ = $config->{expand_summ};

my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );

foreach my $categ (@categs) {
  # Get random selection of sites with predicted rates
  # Must modify to include category
  print "Subsetting $categ sites...\n";
  my $samplecmd = "awk -v categ=\"$categ\" 'BEGIN {srand()} !/^\$/ { if (rand() <= .005) print \$0\"\\t\"0\"\\t\"categ}' $parentdir/output/predicted/chr*.${categ}.txt > $parentdir/output/predicted/${categ}.sub_new.txt";
  &forkExecWait($samplecmd);
}

print "Combining and sorting data...\n";
# Combine null sites with DNMs
my $sortcmd = "cat $parentdir/output/predicted/*.sub_new.txt $parentdir/reference_data/DNMs/*.anno.txt | sort -V -k1,2 > $parentdir/output/rocdat.sort_new.txt";
&forkExecWait($sortcmd);

# Annotate with rate tables (can omit and do this in R after selecting 1M sites)
# for faster processing
# Also need to specify how {group}_7bp_rates.txt tables are generated
# CMD to annotate with common rates
# perl data_mgmt/process_dnms/anno_rate.pl --adj 3 --in /net/bipolar/jedidiah/mutation/output/rocdat.sort.txt --rates /net/bipolar/jedidiah/mutation/common_7bp_rates.txt --seq --out /net/bipolar/jedidiah/mutation/output/rocdat.7bp.1.txt

# Annotate with ERV rates
# perl data_mgmt/process_dnms/anno_rate.pl --adj 3 --in /net/bipolar/jedidiah/mutation/output/rocdat.7bp.1.txt --rates /net/bipolar/jedidiah/mutation/ERV_7bp_rates.txt --out /net/bipolar/jedidiah/mutation/output/rocdat.7bp.txt

# Annotate with AV rates
# perl data_mgmt/process_dnms/anno_rate.pl --adj 3 --in /net/bipolar/jedidiah/mutation/output/rocdat.7bp.1.txt --rates /net/bipolar/jedidiah/mutation/av_7bp_rates.txt --out /net/bipolar/jedidiah/mutation/output/rocdat.7bp.2.txt

##############################################################################
# fork-exec-wait subroutine
##############################################################################
sub forkExecWait {
  my $cmd = shift;
  #print "forkExecWait(): $cmd\n";
  my $kidpid;
  if ( !defined($kidpid = fork()) ) {
    die "Cannot fork: $!";
  }
  elsif ( $kidpid == 0 ) {
    exec($cmd);
    die "Cannot exec $cmd: $!";
  }
  else {
    waitpid($kidpid,0);
  }
}
