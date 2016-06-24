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

# Set options and inputs
my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";

for(i in 1:5000000){
  grep -w "9996" *.dp
}
