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
my $allfiles="$parentdir/output/glf_depth/glf_filelist.txt";

my $chrfiles="$parentdir/output/glf_depth/chr${chr}_glf_filelist.txt";
# my $chrfiles;

my $chrFH;
open($chrFH, '>', $chrfiles) or
  die "Unable to open file $chrfiles : $!";
close($chrFH) or die "Unable to close file: $chrfiles $!";
my $getchrfiles=`grep -w \"chr$chr\" $allfiles > $chrfiles`;

my $chrfilesub;
if($subset==0){
  $numind=3765;
  $chrfilesub=$chrfiles;
} else {

  # Subset file list to 10% of samples
  $chrfilesub="$parentdir/output/glf_depth/chr${chr}_glf_filelist.sub.txt";

  my $samples="/net/bipolar/lockeae/final_freeze/list.txt";
  my $subsetsamples="$parentdir/output/glf_depth/list_sub.txt";
  my $sscmd="cat $samples | perl -ne 'print \$_ if 0.1 > rand;' > $subsetsamples";
  print "Subsetting samples: $sscmd\n";
  &forkExecWait($sscmd);
  my $sscmd2="cat $chrfiles | grep -Fwf $subsetsamples > $chrfilesub";
  print "Subsetting file list: $sscmd2\n";
  &forkExecWait($sscmd2);
  $numind=`wc -l $subsetsamples | cut -d" " -f1`;
  chomp($numind);
  # my $subsetcmd="cat $chrfilesfull | perl -ne 'print $_ if 0.1 > rand;' > $chrfiles";
  # $numind=``
}

# Count total records and specify number of jobs
my $numrecords = `wc -l $chrfilesub | cut -d" " -f1`;
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
