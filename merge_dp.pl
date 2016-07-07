#!/usr/local/bin/perl

##############################################################################
# Adds depth info to motif file
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

my $chr=22;

my $parentdir="/net/bipolar/jedidiah/mutation";
my $dpdir="$parentdir/output/glf_depth/meandp";

my @files = glob("$dpdir/chr$chr.*.txt");

print join("\n", @files);
