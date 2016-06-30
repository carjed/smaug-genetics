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
my $chr=22;
# my $numind=40;
my $subset=0;
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

# initialize and run sbatch file
my $workerbatch = "$parentdir/smaug-genetics/data_mgmt/slurm_process_glfs.$chr.txt";
open my $wFH, '>', $workerbatch or die "can't write to $workerbatch: $!\n";
print $wFH "#!/bin/sh \n";
print $wFH "#SBATCH --mail-type=FAIL \n";
print $wFH "#SBATCH --mail-user=jedidiah\@umich.edu \n";
print $wFH "#SBATCH --ntasks=1 \n";
print $wFH "#SBATCH --mem=2000 \n";
print $wFH "#SBATCH --time 20:00:00 \n";
print $wFH "#SBATCH --job-name=chr${chr}_process_glfs \n";
print $wFH "#SBATCH --partition=nomosix \n";
print $wFH "#SBATCH --array=1-$numjobs \n"; # change to 1-$numjobs
print $wFH "#SBATCH --output=\"$parentdir/output/slurm/slurmJob-%J.out\" --error=\"$parentdir/output/slurm/slurmJob-%J.err\" \n";
print $wFH "srun perl $parentdir/smaug-genetics/data_mgmt/process_glf_worker.pl --chr $chr --ind \${SLURM_ARRAY_TASK_ID} --chunk $chunksize --filelist $chrfilesub\n";
close($wFH) or die "Unable to close file: $workerbatch $!";

my $slurmcmd="sbatch $workerbatch";
# &forkExecWait($slurmcmd);

# my $jobIDfile="$parentdir/output/glf_depth/chr$chr.jobID";
my $rawID=`squeue -u jedidiah | awk 'NR>1 {print \$1}'`;
my $ID=substr($rawID, 0, index($rawID, '_'));
my $datestring = gmtime();
print "Batch job $ID queued at $datestring...\n";

# Continuous loop checks for batch job to complete before proceeding
$datestring = gmtime();
print "Validation started at $datestring...\n";

my $cflag=0;
OUTER:
while($cflag!=1){

  my $logfile="$parentdir/output/glf_depth/chr$chr.data.log2";
  print "Reading $logfile\n";
  # my $logcmd="sacct -j $ID --format=JobID,State | awk 'NR>2 {print \$2}' | sort | uniq | paste -d- -s > $logfile";
  # print "Running $logcmd\n";
  # &forkExecWait($logcmd);

  open my $logFH, '<', $logfile or die "can't open $logfile: $!";

  foreach my $line (<$logFH>){
    chomp($line);
    print "$line\n";
    print "$line\n";
    print "$cflag\n";
    if($line =~ /^COMPLETE/){
      $cflag=1;
      print "VALIDATION COMPLETE\n";
      last;
    }
  }

  close($logFH) or die "Unable to close file: $logfile $!";
  sleep 30;
}

my %filehash=();
my $f_dirlist = "$parentdir/output/glf_depth/chr${chr}_glf_dirlist.txt";
my $getdirlist = "find $parentdir/output/glf_depth/chr$chr -mindepth 1 -maxdepth 1 -type d > $f_dirlist";

my $numdirs = `wc -l $f_dirlist | cut -d" " -f1`;
chomp($numdirs);

# initialize and run sbatch file
my $meandpbatch = "$parentdir/smaug-genetics/data_mgmt/slurm_glf_meandp.$chr.txt";
open my $mdFH, '>', $meandpbatch or die "can't write to $meandpbatch: $!\n";
print $mdFH "#!/bin/sh \n";
print $mdFH "#SBATCH --mail-type=FAIL \n";
print $mdFH "#SBATCH --mail-user=jedidiah\@umich.edu \n";
print $mdFH "#SBATCH --ntasks=1 \n";
print $mdFH "#SBATCH --mem=2000 \n";
print $mdFH "#SBATCH --time 20:00:00 \n";
print $mdFH "#SBATCH --job-name=chr${chr}_glf_meandp \n";
print $mdFH "#SBATCH --partition=nomosix \n";
print $mdFH "#SBATCH --array=1-$numdirs \n"; # change to 1-$numjobs
print $mdFH "#SBATCH --output=\"$parentdir/output/slurm/slurmJob-%J.out\" --error=\"$parentdir/output/slurm/slurmJob-%J.err\" \n";
print $mdFH "srun perl $parentdir/smaug-genetics/data_mgmt/process_glf_meandp.pl --chr $chr --ind \${SLURM_ARRAY_TASK_ID}";
close($mdFH) or die "Unable to close file: $meandpbatch $!";

$slurmcmd="sbatch $meandpbatch";
&forkExecWait($slurmcmd);

$rawID=`squeue -u jedidiah | awk 'NR>1 {print \$1}'`;
$ID=substr($rawID, 0, index($rawID, '_'));
$datestring = gmtime();
print "Batch job $ID queued at $datestring...\n";

sub is_folder_empty {
    my $dirname = shift;
    opendir(my $dh, $dirname) or die "Not a directory";
    return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
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
