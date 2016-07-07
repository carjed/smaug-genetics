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

my @chrs=(16..18, 4);
foreach my $i (@chrs){
  print "Processing chr$i...\n";
  my $pgmcmd="perl /net/bipolar/jedidiah/mutation/smaug-genetics/data_mgmt/process_glf_master.pl --chr $i";
  &forkExecWait($pgmcmd);
  print "==========================\n";
}

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
