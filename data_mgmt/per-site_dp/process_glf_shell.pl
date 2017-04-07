#!/usr/local/bin/perl

##############################################################################
# Master script for glf processing
##############################################################################

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Math::Round;
use Cwd;
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

my @chrs=(17..18, 4);
foreach my $i (@chrs){
  print "Processing chr$i...\n";
  my $pgmcmd="perl $parentdir/smaug-genetics/data_mgmt/per-site_dp/process_glf_master.pl --chr $i";
  forkExecWait($pgmcmd);
  print "==========================\n";
}
