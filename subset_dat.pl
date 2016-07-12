#!/usr/local/bin/perl

##############################################################################
# Subset data for logit model
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

my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";
my $categ="AT_CG";

# Index motif file names
my $f_mlist = "$parentdir/output/7bp_1000k_rates.txt";
open my $mlist, '<', $f_mlist or die "can't open $f_mlist: $!";

our %fhash=();
my @fn;
while(<$mlist>){
  chomp;
  my @line=split(/\t/, $_);
  my $seq=$line[1];
  my $cat=$line[2];

  if($cat eq $categ){
    #my $key=join("\t", @line[0 .. 1]);
    #my $pcs=join("\t", @line[2 .. $#line]);

    my $filename="$parentdir/output/logmod_data/${categ}_tmp_$seq.txt";
    push(@fn, $filename);
    $fhash{$seq}=$filename;
    # print "$hash{$_}\n";
  }
}

my %handles = get_write_handles(@fn);

# foreach(@fn){
#   my $ftest=$handles{$_};
#   print "$ftest\n";
# }
#
# foreach(values %hash){
#   print "$_\n";
#   print {$handles{$_}} "$_\n";
# }

# Read full data and subset by motif
my $f_positions = "$parentdir/output/logmod_data/${categ}_full.txt";
open my $positions, '<', $f_positions or die "can't open $f_positions: $!";

while(<$positions>){
  chomp;
  my @line=split(/\t/, $_);
  my $motif=$line[3];
  # print "$motif\n";

  # if(exists $hash{$motif}){
    my $file=$fhash{$motif};
    print {$handles{$file}} "$_\n";
  # }
}

# Subroutine reads array of filenames and returns file handles
sub get_write_handles {
  my @file_names = @_;
  my %file_handles;
  foreach (@file_names) {
    open my $fh, '>>', $_ or next;
    $file_handles{$_} = $fh;
  }
  return %file_handles;
}
