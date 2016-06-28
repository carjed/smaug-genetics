#!/usr/local/bin/perl

##############################################################################
# Script scans glf file and extracts only the relevant info (chr, pos, ref, dp)
# for every 10th base
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

my $parentdir="/net/bipolar/jedidiah/mutation";

my $chr=3;

my $filelist="$parentdir/output/glf_depth/glf_filelist.txt";
my $chrfilelist="$parentdir/output/glf_depth/chr${chr}_glf_filelist.txt";

my $getchrfiles="grep -w \"chr$chr\"  $filelist > $chrfilelist";
my $count = `wc -l < $chrfilelist`;
print "$count\n";

my $jobIDfile="$parentdir/output/glf_depth/chr$chr.jobID";
my $getjobID=`squeue -u jedidiah | awk '{print \$1}'`;
print "$getjobID\n";
# &forkExecWait($getjobID);

# open my $jobID, '<', $jobIDfile or die "can't open $jobIDfile: $!";
# my $ID;
# while(<$jobID>){
#   chomp;
#   $ID=$_;
# }
#
# my $logfile="$parentdir/output/glf_depth/chr$chr.data.log";
# my $logcmd="sacct -j $ID --format=JobID,State | awk 'NR>2 {print \$2}' | sort | uniq -c > $logfile";
# &forkExecWait($logcmd);
#
# open my $log, '<', $logfile or die "can't open $logfile: $!";
# while(<$log>){
#   if $_
# }

##############################################################################
# fork-exec-wait subroutine
##############################################################################
sub forkExecWait {
    my $cmd = shift;
    #print "forkExecWait(): $cmd\n";
    my $kidpid;
    if ( !defined($kidpid = fork()) ) {
		    die "Cannot fork: $!";
    }
    elsif ( $kidpid == 0 ) {
		    exec($cmd);
		    die "Cannot exec $cmd: $!";
    }
    else {
		    waitpid($kidpid,0);
    }
}
