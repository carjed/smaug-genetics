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

my $chr=4;
my $numsamples=400;
my $parentdir="/net/bipolar/jedidiah/mutation";
  make_path("$parentdir/output/glf_depth/chr$chr");

my $filelist="$parentdir/output/glf_depth/glf_filelist.txt";
my $chrfilelist="$parentdir/output/glf_depth/chr${chr}_glf_filelist.txt";
my $chrfile;
open($chrfile, '>', $chrfilelist) or
die "Unable to open file $chrfilelist : $!";
close($chrfile) or die "Unable to close file: $chrfilelist $!";
my $getchrfiles=`grep -w \"chr$chr\" $filelist > $chrfilelist`;

my $count = `wc -l $chrfilelist | cut -d" " -f1`;
die "wc failed: $?" if $?;
chomp($count);
my $njobs=ceil($count/400);
print "Number of records to process in chr$chr: $count\n";

# initialize output
my $outfile = "$parentdir/smaug-genetics/data_mgmt/slurm_process_glfs.$chr.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";
print OUT "#!/bin/sh \n";
print OUT "#SBATCH --mail-type=ALL \n";
print OUT "#SBATCH --mail-user=jedidiah\@umich.edu \n";
print OUT "#SBATCH --ntasks=1 \n";
print OUT "#SBATCH --mem=10000 \n";
print OUT "#SBATCH --time 10:00:00 \n";
print OUT "#SBATCH --job-name=chr${chr}_process_glfs \n";
print OUT "#SBATCH --partition=nomosix \n";
print OUT "#SBATCH --array=1-1 \n";
print OUT "#SBATCH --output=\"$parentdir/output/slurm/slurmJob-%J.out\" --error=\"$parentdir/output/slurm/slurmJob-%J.err\" \n";
print OUT "srun perl $parentdir/smaug-genetics/data_mgmt/process_glf.pl --chr $chr --ind \${SLURM_ARRAY_TASK_ID}\n";
close(OUT) or die "Unable to close file: $outfile $!";

my $slurmcmd="sbatch $outfile";
&forkExecWait($slurmcmd);

my $jobIDfile="$parentdir/output/glf_depth/chr$chr.jobID";
my $rawID=`squeue -u jedidiah | awk 'NR>1 {print \$1}'`;
my $ID=substr($rawID, 0, index($rawID, '_'));
print "$ID\n";
# &forkExecWait($getjobID);

my $f_dirlist = "$parentdir/output/glf_depth/chr${chr}_glf_dirlist.txt";
my $getdirlist = "find $parentdir/output/glf_depth/chr$chr -mindepth 1 -maxdepth 1 -type d > $f_dirlist";
&forkExecWait($getdirlist);
open my $dirlist, '<', $f_dirlist or die "can't open $f_dirlist: $!";

while(<$dirlist>){
  print "Validating files in $_...";
  while (1) {
    my $numfiles=`$_/*.dp | wc -l`;
    my $chknum=`$_/*.ok | wc -l`;
    last if $chknum==$numsamples;
    sleep 1;
  }
  # run glf_depth.pl
  print "COMPLETE\n";
}

# my $logfile="$parentdir/output/glf_depth/chr$chr.data.log";
# my $logcmd="sacct -j $ID --format=JobID,State | awk 'NR>2 {print \$2}' | sort | uniq -c > $logfile";
# &forkExecWait($logcmd);

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
