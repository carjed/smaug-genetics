#!/usr/local/bin/perl

##############################################################################
# Master script for glf processing
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
my $glfdir = $config->{glfdir};
my $samples = $config->{samples};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

my $filelistcmd = "find $glfdir/chr* -mindepth 1 -type f -name \"*.glf\" | grep -e \"chr[0-9]\" | grep -Fwf $samples > $parentdir/output/glf_depth/glf_filelist.txt";

my @chrs=(1..22);
foreach my $i (@chrs){
  print "Processing chr$i...\n";
  my $pgmcmd="perl $parentdir/smaug-genetics/data_mgmt/per-site_dp/process_glf_master.pl $i";
  forkExecWait($pgmcmd);
  print "==========================\n";
}
