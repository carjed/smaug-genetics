#!/usr/local/bin/perl

##############################################################################
# Master script for glf processing
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
use Benchmark;
use Tie::File;

# Specify parameters
my $chr=4;
# my $numind=40;
my $subset=1;
my $numind=3765;
my $chunksize=40; # No. of records to process in each worker job
my $parentdir="/net/bipolar/jedidiah/mutation";
  make_path("$parentdir/output/glf_depth/chr$chr");
my $filelist="$parentdir/output/glf_depth/glf_filelist.txt";

my $chrfilesfull="$parentdir/output/glf_depth/chr${chr}_glf_filelist.txt";
my $chrfiles;

if($subset==0){
  $numind=3765;
  $chrfiles=$chrfilesfull;
} else {

  # Subset file list to 10% of samples
  $chrfiles="$parentdir/output/glf_depth/chr${chr}_glf_filelist.sub.txt";

  my $samples="/net/bipolar/lockeae/final_freeze/list.txt";
  my $subsetsamples="$parentdir/output/glf_depth/list_sub.txt";
  my $sscmd="cat $samples | perl -ne 'print \$_ if 0.1 > rand;' > $subsetsamples";
  &forkExecWait($sscmd);
  my $sscmd2="cat $chrfilesfull | grep -Fwf $subsetsamples > $chrfiles";
  &forkExecWait($sscmd2);
  $numind=`wc -l $subsetsamples | cut -d" " -f1`;
  # my $subsetcmd="cat $chrfilesfull | perl -ne 'print $_ if 0.1 > rand;' > $chrfiles";
  # $numind=``
}

my $chrfile;
open($chrfile, '>', $chrfiles) or
  die "Unable to open file $chrfiles : $!";
close($chrfile) or die "Unable to close file: $chrfiles $!";
my $getchrfiles=`grep -w \"chr$chr\" $filelist > $chrfiles`;

# Count total records and specify number of jobs
my $numrecords = `wc -l $chrfiles | cut -d" " -f1`;
die "wc failed: $?" if $?;
chomp($numrecords);
my $numjobs=ceil($numrecords/$chunksize);
print "Number of records to process in chr$chr: $numrecords\n";
print "Number of individuals to be processed: $numind\n";

##############################################################################
# fork-exec-wait subroutine
##############################################################################
sub forkExecWait {
  my $cmd = shift;
  #print "forkExecWait(): $cmd\n";
  my $kidpid;
  if ( !defined($kidpid = fork()) ) {
	  die "Cannot fork: $!";
  } elsif ($kidpid==0) {
	  exec($cmd);
	  die "Cannot exec $cmd: $!";
  } else {
	  waitpid($kidpid,0);
  }
}
