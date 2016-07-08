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

# Specify parameters
my $i;
my $categ;
# my $out;
my $help=0;
my $man=0;

my $parentdir="/net/bipolar/jedidiah/mutation";
my $dpdir="$parentdir/output/glf_depth/meandp";

# Set options and inputs
GetOptions ('in=i'=> \$i,
# 'out=s'=> \$out,
'categ=s' => \$categ,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;

my $in=`ls -d $parentdir/output/logmod_data/motifs/$categ/* | sed -n '$i'p`;
chomp($in);
# my $new = $old =~ s/foo/bar/r;
my $out = $in =~ s/.txt/_dp.txt/r;
$out =~ s/$categ\//$categ\/dp\//;

my $inFH;
open($inFH, '<', $in) or
  die "Unable to open file $in : $!";

my $curchr=0;
my %dphash=();

# my $ = "$parentdir/smaug-genetics/data_mgmt/slurm_process_glfs.$chr.txt";
open my $outFH, '>', $out or die "can't write to $out: $!\n";

while(<$inFH>){
  chomp;
  my @line=split(/\t/, $_);
  my $pos=$line[0];
  my $chr=$line[1];

  if($chr!=$curchr){
    $curchr=$chr;
    %dphash=();

    my $chrfile="$dpdir/chr$curchr.dp";
    my $chrFH;
    open($chrFH, '<', $chrfile) or
      die "Unable to open file $chrfile : $!";

    while(my $dp=<$chrFH>){
      chomp($dp);
      my @dpline=split(/\t/, $dp);
      my $dp_pos=$dpline[0];
      my $depth=$dpline[1];

      $dphash{$dp_pos}=$depth;
    }
  }

  my $poslim=rounddown($pos,10);

  if(exists($dphash{$poslim})){
    my $dpout=$dphash{$poslim};
    print $outFH "$_\t$dpout\n";
  }
}


##############################################################################
# Rounding subroutines
##############################################################################
sub roundup
{
  my $num = shift;
  my $roundto = shift || 1;

  return int(ceil($num/$roundto))*$roundto;
}

sub rounddown
{
  my $num = shift;
  my $roundto = shift || 1;

  return int(floor($num/$roundto))*$roundto;
}
