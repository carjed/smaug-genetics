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

my $numind=40;
# my $numind=3765;

# Number of records to process in each worker job
my $numsamples=40;

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
my $njobs=ceil($count/$numsamples);
print "Number of records to process in chr$chr: $count\n";

# initialize output
my $outfile = "$parentdir/smaug-genetics/data_mgmt/slurm_process_glfs.$chr.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";
print OUT "#!/bin/sh \n";
print OUT "#SBATCH --mail-type=FAIL \n";
print OUT "#SBATCH --mail-user=jedidiah\@umich.edu \n";
print OUT "#SBATCH --ntasks=1 \n";
print OUT "#SBATCH --mem=10000 \n";
print OUT "#SBATCH --time 10:00:00 \n";
print OUT "#SBATCH --job-name=chr${chr}_process_glfs \n";
print OUT "#SBATCH --partition=nomosix \n";
print OUT "#SBATCH --array=1-1 \n";
print OUT "#SBATCH --output=\"$parentdir/output/slurm/slurmJob-%J.out\" --error=\"$parentdir/output/slurm/slurmJob-%J.err\" \n";
print OUT "srun perl $parentdir/smaug-genetics/data_mgmt/process_glf.pl --chr $chr --ind \${SLURM_ARRAY_TASK_ID} --chunk 40 \n";
close(OUT) or die "Unable to close file: $outfile $!";

my $slurmcmd="sbatch $outfile";
&forkExecWait($slurmcmd);

my $jobIDfile="$parentdir/output/glf_depth/chr$chr.jobID";
my $rawID=`squeue -u jedidiah | awk 'NR>1 {print \$1}'`;
my $ID=substr($rawID, 0, index($rawID, '_'));
my $datestring = gmtime();
print "Batch job $ID started at $datestring...\n";

my %filehash=();
my $f_dirlist = "$parentdir/output/glf_depth/chr${chr}_glf_dirlist.txt";
my $getdirlist = "find $parentdir/output/glf_depth/chr$chr -mindepth 1 -maxdepth 1 -type d > $f_dirlist";

$datestring = gmtime();
print "Validation started at $datestring...\n";
my $cflag=0;
while($cflag!=1){

  # To-do:
  # [X] move this code chunk inside first while loop above
  # [X] run glf_depth.pl
  # [-] delete glfs and OK file in subdir once mean depth file finished
  # [X] mark completed directories; skip when re-scanning (create hash table!)
  my $getdirlist = "find $parentdir/output/glf_depth/chr$chr -mindepth 1 -maxdepth 1 -type d > $f_dirlist";
  &forkExecWait($getdirlist);
  open my $dirlist, '<', $f_dirlist or die "can't open $f_dirlist: $!";

  while(<$dirlist>){
    chomp;
    if(exists($filehash{$_}) && $filehash{$_}!=1){
      # print "Validating files in $_...";
      # while (1) {
        # my $numfiles=`ls $_/\*.dp | wc -l | cut -d" " -f1`;
        my @filecount=glob("$_/*.dp");
        my $numfiles = scalar @filecount;
        # chomp($numfiles);
        my $chknum=`wc -l $_/samples.ok | cut -d" " -f1`;
        chomp($chknum);
        if($chknum==$numind){
          # run glf_depth.pl on chunk

          my $perlcmd = "perl $parentdir/smaug-genetics/data_mgmt/glf_depth.pl --chr $chr --dir $_";
          &forkExecWait($perlcmd);
          $filehash{$_}=1;

          my $meandp = "$parentdir/output/glf_depth/chr$chr.1.5000000.txt"
          if(-e ){
            my $rmcmd="rm -f $_/\*.dp"
          }
          # print "COMPLETE\n";
          # last;
        }
        # sleep 1;
      # }

    } else {
      $filehash{$_}=0;
    }
  }

  my $logfile="$parentdir/output/glf_depth/chr$chr.data.log";
  my $logcmd="sacct -j $ID --format=JobID,State | awk 'NR>2 {print \$2}' | sort | uniq | paste -d- -s > $logfile";
  &forkExecWait($logcmd);
  open my $log, '<', $logfile or die "can't open $logfile: $!";

  # To-do: also validate that mean depth files created for each subdir before script finishes
  while(<$log>){
    chomp;
    if($_ eq "COMPLETED"){
      $cflag=1;
      print "VALIDATION COMPLETE\n";
      last;
    }
  }

  sleep 300;
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
    }
    elsif ( $kidpid == 0 ) {
		    exec($cmd);
		    die "Cannot exec $cmd: $!";
    }
    else {
		    waitpid($kidpid,0);
    }
}
