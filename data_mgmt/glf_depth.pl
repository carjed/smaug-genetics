#!/usr/local/bin/perl

##############################################################################
#
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
my $chr=1;

my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";
my $glfdir="$parentdir/output/glf_depth/chr$chr/";

print "$glfdir\n";

opendir (DIR, $glfdir) or die $!;

my @dirs = grep {-d "$glfdir/$_" && ! /^\.{1,2}$/} readdir(DIR);
print "$_\n" foreach @dirs;

my @range = split(/./, @dirs[0]);
my $start = $range[0];
my $end = $range[1];

print "$start\t$end\n";

# for(i in 1:5000000){
#   grep -w "9996" *.dp
# }

my %hash=();
my %hashn=();

# $hashn{$_}=0 for

my $outfile = "$parentdir/output/glf_depth/test.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

# glob ('/path/to/dir/*');
my @files = glob("$parentdir/output/glf_depth/1497-RMM-401*");
foreach my $file (@files) {
  print $file . "\n";
  open my $sample, '<', $file or die "can't open $file: $!";

  while (<$sample>){
  	chomp;
  	my @line=split(/\t/, $_);
  	my $pos=$line[1];
  	# my $vals=join("\t", nearest(0.0001, @line[1 .. $#line]));
  	my $dp=$line[3];
  	$hash{$pos}+=$dp;
    $hashn{$pos}+=1;
  }

}

$hash{$_}=$hash{$_}/$hashn{$_} foreach (keys%hash);

print OUT "$_\t$hash{$_}\n" foreach (sort {$a <=> $b} (keys%hash));


sub roundup
{
  my $num = shift;
  my $roundto = shift || 1;

  return int(ceil($num/$roundto))*$roundto;
}
