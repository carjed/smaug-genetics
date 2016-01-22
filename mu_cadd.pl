#!/usr/local/bin/perl

##############################################################################
# Merge predicted mutation rates with C-scores
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
use Compress::Zlib;
use IO::Compress::Gzip;

# Set options and inputs
my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";

my $f_mus = "/net/bipolar/jedidiah/mutation/output/predicted/AT_GC_out.txt"; #main line for full processing
open my $mus, '<', $f_mus or die "can't open $f_mus: $!";

my $f_cadd = "/net/dumbo/home/lockeae/CADD/whole_genome_SNVs.tsv.gz";
open my $cadd, '<', $f_cadd or die "can't open $f_cadd: $!";

my $vcf = gzopen($f_cadd, "rb") or
  die "can't open $invcf: $gzerrno";

while(<$mus>){
  chomp;
  my @line=split(/\t/, $_);
  my $CHR=$line[0];
  my $POS=$line[1];
  my $MU=$line[2];
}
