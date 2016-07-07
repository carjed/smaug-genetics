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

my $chr=21;

my $parentdir="/net/bipolar/jedidiah/mutation";
my $dpdir="$parentdir/output/glf_depth/meandp";

# my $rawfiles = `ls -v $dpdir/chr$chr.*.txt | xargs cat >> $dpdir/chr$chr.dp`;
my $rawfiles = `ls -v $dpdir/chr$chr.*.txt`;

# print "$rawfiles\n";
my @files = split(/\n/, $rawfiles);

my $out="$dpdir/chr$chr.dp";
open my $outFH, '>', $out or die "can't write to $out: $!\n";
while(<@files>){
  my $blockFH;
  open($blockFH, '>', $_) or
    die "Unable to open file $_ : $!";

  my @filepath=split m%/%, $_;

  # print join("\n", @filepath);
  my @range=split(/\./, $filepath[-1]);
  # print join("\n", @range);

  my $start=$range[1];
  my $end=$range[2];

  print "$start\t$end\n";
  while(my $line=<$blockFH>){
    chomp;
    my @elements=split(/\t/, $line);
    my $pos=$elements[0];

    if($pos>=$start && $pos<=$end){
      print $outFH "$line\n";
    }
  }
}


# my @sortfiles=sort by_number @files;
# print join("\n", @files);

# sub by_number {
#     my ( $anum ) = $a =~ /(\d+)/;
#     my ( $bnum ) = $b =~ /(\d+)/;
#     ( $anum || 0 ) <=> ( $bnum || 0 );
# }
