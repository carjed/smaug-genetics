#!/usr/local/bin/perl

##############################################################################
# Script scans glf file and extracts only the relevant info (chr, pos, ref, dp)
# for every 10th base
##############################################################################
use strict;
use warnings;
use POSIX;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

my $parentdir = $config->{parentdir};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

# Set options and inputs
my $index;
my $chr;
my $chunksize=400;
my $filelist;

GetOptions ('ind=i'=> \$index,
'chr=i' => \$chr,
'chunk=i' => \$chunksize,
'filelist=s' => \$filelist);

# my $filelist="$parentdir/output/glf_depth/chr${chr}_glf_filelist.txt";
open my $files, '<', $filelist or die "$filelist: $!";
my $NUMFILES=2217585;

print "Worker ID: $index\n";

my $start=($index-1)*$chunksize+1;
my $end=$index*$chunksize;

my @filerange;
while(<$files>) {
  chomp;
  if (($. == $start) .. ($. == $end)) {
    push @filerange, $_;
  }
}

my $numruns=scalar @filerange;
# print "Running loop on $numruns\n";
foreach my $sample (@filerange){
  my @filepath=split m%/%, $sample;
  my $froot="$filepath[8]/$filepath[9]";
  my $fname="$froot/$filepath[10].dp";
  make_path("$parentdir/output/glf_depth/$froot");

  my $file="$parentdir/output/glf_depth/$fname";
  if(-e $file){
    my $skip=1;
  } else {
    my $glfcmd="samtools-hybrid glfview $sample | cut -f1-4 | awk '\$2%10==0 && \$3 ~ /[ACGT]/' > $parentdir/output/glf_depth/$fname";
    # print "$glfcmd\n";
    forkExecWait($glfcmd);
  }

  my $okfile="$parentdir/output/glf_depth/$froot/samples.ok";
  open(OUT, '>>', $okfile) or die "can't write to $okfile: $!\n";
  print OUT "$fname: OK\n";
  close(OUT) or die "Unable to close file: $okfile $!";
  # print "$sample\n";
}
