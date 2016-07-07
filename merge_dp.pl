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

my $rawfiles = `ls -v $dpdir/chr$chr.*.txt`;

my @files = split($rawfiles, "\n");



# my @sortfiles=sort by_number @files;
print join("\n", @files);

# sub by_number {
#     my ( $anum ) = $a =~ /(\d+)/;
#     my ( $bnum ) = $b =~ /(\d+)/;
#     ( $anum || 0 ) <=> ( $bnum || 0 );
# }
