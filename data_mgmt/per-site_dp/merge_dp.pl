#!/usr/local/bin/perl

##############################################################################
# Merge chr block files
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

# my $chr=21;
# my $help=0;
# my $man=0;
#
# # Set options and inputs
# GetOptions ('chr=i'=> \$chr,
# 'help|?'=> \$help,
# man => \$man) or pod2usage(1);

my $parentdir="/net/bipolar/jedidiah/mutation";
my $dpdir="$parentdir/output/glf_depth/meandp";
# my @chrs=(1..3, 5..17, 19, 20);
my @chrs=(18);
foreach my $chr (@chrs){
  # my $rawfiles = `ls -v $dpdir/chr$chr.*.txt | xargs cat >> $dpdir/chr$chr.dp`;
  my $rawfiles = `ls -v $dpdir/chr$chr.*.txt`;

  # print "$rawfiles\n";
  my @files = split(/\n/, $rawfiles);

  my $out="$dpdir/chr$chr.dp";
  open my $outFH, '>', $out or die "can't write to $out: $!\n";
  while(<@files>){
    my $blockFH;
    open($blockFH, '<', $_) or
      die "Unable to open file $_ : $!";

    my @filepath=split m%/%, $_;

    # print join("\n", @filepath);
    my @range=split(/\./, $filepath[-1]);
    # print join("\n", @range);

    my $start=$range[1];
    my $end=$range[2];

    print "$start\t$end\n";
    while(my $line=<$blockFH>){
      chomp($line);
      my @elements=split(/\t/, $line);
      my $pos=$elements[0];

      if($pos>=$start && $pos<=$end){
        print $outFH "$line\n";
      }
    }
  }
}
