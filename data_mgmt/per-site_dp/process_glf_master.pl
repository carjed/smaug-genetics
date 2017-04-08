#!/usr/local/bin/perl

##############################################################################
# Master script for glf processing
##############################################################################
use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Path qw(make_path);
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

my $email = $config->{email};
my $parentdir = $config->{parentdir};
my $samples = $config->{samples};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

# Get chromosome from args
my $chr=$ARGV[0];

# my $numind=40;
my $subset=0;
my $numind=3765;
my $chunksize=40; # No. of records to process in each worker job

make_path("$parentdir/output/glf_depth/chr$chr");
my $allfiles="$parentdir/output/glf_depth/glf_filelist.txt";
my $chrfiles="$parentdir/output/glf_depth/chr${chr}_glf_filelist.txt";

my $today = POSIX::strftime('%Y%m%d', localtime);
my $slurmdir = "$parentdir/output/slurm/$today";
  make_path("$slurmdir");

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
  # Not used--tested if we can maintain accuracy and reduce workload
  # by subsampling file list to 10% of samples
  $chrfilesub="$parentdir/output/glf_depth/chr${chr}_glf_filelist.sub.txt";

  my $subsetsamples="$parentdir/output/glf_depth/list_sub.txt";
  my $sscmd="cat $samples | perl -ne 'print \$_ if 0.1 > rand;' > $subsetsamples";
  print "Subsetting samples: $sscmd\n";
  forkExecWait($sscmd);
  my $sscmd2="cat $chrfiles | grep -Fwf $subsetsamples > $chrfilesub";
  print "Subsetting file list: $sscmd2\n";
  forkExecWait($sscmd2);
  $numind=`wc -l $subsetsamples | cut -d" " -f1`;
  chomp($numind);
}

# Count total records and specify number of jobs
my $numrecords = `wc -l $chrfilesub | cut -d" " -f1`;
chomp($numrecords);

my $numjobs=ceil($numrecords/$chunksize);
print "Number of records to process in chr$chr: $numrecords\n";
print "Number of individuals to be processed: $numind\n";

# initialize and run sbatch file
my $jobcmd="chr${chr}_process_glfs";
my $workerbatch = "$parentdir/slurm/slurm_process_glfs.$chr.txt";
open my $wFH, '>', $workerbatch or die "can't write to $workerbatch: $!\n";
print $wFH "#!/bin/sh \n";
print $wFH "#SBATCH --mail-type=FAIL \n";
print $wFH "#SBATCH --mail-user=$email \n";
print $wFH "#SBATCH --ntasks=1 \n";
print $wFH "#SBATCH --mem=2000 \n";
print $wFH "#SBATCH --time 20:00:00 \n";
print $wFH "#SBATCH --job-name=$jobcmd \n";
print $wFH "#SBATCH --partition=nomosix \n";
print $wFH "#SBATCH --array=1-$numjobs \n"; # change to 1-$numjobs
print $wFH "#SBATCH --requeue \n";
print $wFH "#SBATCH --exclude=r30,r14,inpsyght \n";
# print $wFH "#SBATCH --exclude=psoriasis-mc01,psoriasis-mc02 \n";
print $wFH "#SBATCH --output=\"$slurmdir/slurmJob-%J.out\" --error=\"$slurmdir/slurmJob-%J.err\" \n";
print $wFH "srun perl $parentdir/smaug-genetics/data_mgmt/per-site_dp/process_glf_worker.pl --chr $chr --ind \${SLURM_ARRAY_TASK_ID} --chunk $chunksize --filelist $chrfilesub\n";
close($wFH) or die "Unable to close file: $workerbatch $!";

my $slurmcmd="sbatch $workerbatch";
forkExecWait($slurmcmd);

# my $jobIDfile="$parentdir/output/glf_depth/chr$chr.jobID";
my $rawID=`squeue --format \"%.24i %.9P %.24j %.8u %.2t %.10M %.6D\" -n $jobcmd | awk 'NR==2 {print \$1}'`;
my $ID=substr($rawID, 0, index($rawID, '_'));
my $datestring = localtime();
print "Batch job $ID queued at $datestring...\n";

&validate_slurm($ID, $numjobs, 120);

my %filehash=();
my $f_dirlist = "$parentdir/output/glf_depth/chr${chr}_glf_dirlist.txt";
my $getdirlist = "find $parentdir/output/glf_depth/chr$chr -mindepth 1 -maxdepth 1 -type d > $f_dirlist";
forkExecWait($getdirlist);

my $numdirs = `wc -l $f_dirlist | cut -d" " -f1`;
chomp($numdirs);

# initialize and run sbatch file
my $x = 1000 + int(rand(9999 - 1000));
$jobcmd="chr${chr}_glf_meandp_$x";
my $meandpbatch = "$parentdir/slurm/slurm_glf_meandp.$chr.txt";
open my $mdFH, '>', $meandpbatch or die "can't write to $meandpbatch: $!\n";
print $mdFH "#!/bin/sh \n";
print $mdFH "#SBATCH --mail-type=FAIL \n";
print $mdFH "#SBATCH --mail-user=$email \n";
print $mdFH "#SBATCH --ntasks=1 \n";
print $mdFH "#SBATCH --mem=2000 \n";
print $mdFH "#SBATCH --time 20:00:00 \n";
print $mdFH "#SBATCH --job-name=$jobcmd \n";
print $mdFH "#SBATCH --partition=nomosix \n";
print $mdFH "#SBATCH --array=1-$numdirs \n";
print $mdFH "#SBATCH --requeue \n";
print $mdFH "#SBATCH --exclude=r30,r14,inpsyght \n";
# print $mdFH "#SBATCH --exclude=psoriasis-mc01,psoriasis-mc02 \n";
print $mdFH "#SBATCH --output=\"$slurmdir/slurmJob-%J.out\" --error=\"$slurmdir/slurmJob-%J.err\" \n";
print $mdFH "srun perl $parentdir/smaug-genetics/data_mgmt/per-site_dp/process_glf_meandp.pl --chr $chr --ind \${SLURM_ARRAY_TASK_ID}";
close($mdFH) or die "Unable to close file: $meandpbatch $!";

$slurmcmd="sbatch $meandpbatch";
forkExecWait($slurmcmd);

$rawID=`squeue --format \"%.24i %.9P %.24j %.8u %.2t %.10M %.6D\" -n $jobcmd | awk 'NR==2 {print \$1}'`;
$ID=substr($rawID, 0, index($rawID, '_'));
$datestring = localtime();
print "Batch job $ID queued at $datestring...\n";

&validate_slurm($ID, $numdirs, 600);

##############################################################################
# Subroutine checks if folder is empty or not
##############################################################################
sub is_folder_empty {
    my $dirname = shift;
    opendir(my $dh, $dirname) or die "Not a directory";
    return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

##############################################################################
# Subroutine checks jobs in SBATCH array and requeues if failed
##############################################################################
sub validate_slurm {
  my $ID=shift;
  my $numjobs=shift;
  my $interval=shift;

  $datestring = localtime();
  print "Validation started at $datestring...\n";

  sleep 60;
  my %statushash=();

  my $njs=$numjobs-1;
  foreach my $i (1..$njs){
    $statushash{$i}=0;
  }

  my $cflag=0;
  OUTER:
  while($cflag!=1){

    # my $logfile="$parentdir/output/glf_depth/chr$chr.data.log";
    # my $logcmd="sacct -j $ID --format=JobID,State | awk 'NR>2 {print \$2}' | sort | uniq | paste -d- -s > $logfile";
    # forkExecWait($logcmd);
    #
    # open my $logFH, '<', $logfile or die "can't open $logfile: $!";
    my $numcompleted=0;
    foreach my $i (1..$njs){
      # my $jobcmd="perl $parentdir/smaug-genetics/data_mgmt/process_glf_worker.pl --chr $chr --ind $i --chunk $chunksize --filelist $chrfilesub";
      # print $jobFH "$jobcmd\n";

      if($statushash{$i}!=1){
        my $grepstr = "${ID}_$i ";
        my $status=`sacct -X -j $grepstr --format=jobid%30,state | awk 'NR>2 {print \$2}'`;
        chomp($status);

        if($status eq "COMPLETED"){
          $statushash{$i}=1;
        } elsif($status eq "FAILED"){
          $statushash{$i}=0;
          # forkExecWait($jobcmd);
          my $requeuecmd="scontrol requeue $grepstr";
          forkExecWait($requeuecmd);
        } else {
          $statushash{$i}=0;
        }
      }
      $numcompleted = sum values %statushash;

      if($numcompleted>=$njs){
        $cflag=1;
        print "VALIDATION COMPLETE\n";
        last OUTER;
      }
    }

    $datestring = localtime();
    print "$numcompleted of $njs jobs finished at $datestring\n";
    # close($logFH) or die "Unable to close file: $logfile $!";
    sleep $interval;
  }
}

##############################################################################
# Alternative validation scheme--doesn't rely on sacct
##############################################################################
# sub validate_slurm2{
#
# }
