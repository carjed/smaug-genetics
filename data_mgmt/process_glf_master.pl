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

# initialize and run sbatch file
my $outfile = "$parentdir/smaug-genetics/data_mgmt/slurm_process_glfs.$chr.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";
print OUT "#!/bin/sh \n";
print OUT "#SBATCH --mail-type=FAIL \n";
print OUT "#SBATCH --mail-user=jedidiah\@umich.edu \n";
print OUT "#SBATCH --ntasks=1 \n";
print OUT "#SBATCH --mem=2000 \n";
print OUT "#SBATCH --time 20:00:00 \n";
print OUT "#SBATCH --job-name=chr${chr}_process_glfs \n";
print OUT "#SBATCH --partition=nomosix \n";
print OUT "#SBATCH --array=1-$numjobs \n"; # change to 1-$numjobs
print OUT "#SBATCH --output=\"$parentdir/output/slurm/slurmJob-%J.out\" --error=\"$parentdir/output/slurm/slurmJob-%J.err\" \n";
print OUT "srun perl $parentdir/smaug-genetics/data_mgmt/process_glf_worker.pl --chr $chr --ind \${SLURM_ARRAY_TASK_ID} --chunk $chunksize \n";
close(OUT) or die "Unable to close file: $outfile $!";

my $slurmcmd="sbatch $outfile";
&forkExecWait($slurmcmd);

my $jobIDfile="$parentdir/output/glf_depth/chr$chr.jobID";
my $rawID=`squeue -u jedidiah | awk 'NR>1 {print \$1}'`;
my $ID=substr($rawID, 0, index($rawID, '_'));
my $datestring = gmtime();
print "Batch job $ID queued at $datestring...\n";

# Continuous loop checks if all files in each 5Mb subdirectory are present
# If so, runs script to get mean per position, deletes files when done
# To-do:
# [-] delete glfs and OK file in subdir once mean depth file finished
# [-] Fix $chknum==$numind equivalency--currently throws error if slurm job is pending
my %filehash=();
my $f_dirlist = "$parentdir/output/glf_depth/chr${chr}_glf_dirlist.txt";
my $getdirlist = "find $parentdir/output/glf_depth/chr$chr -mindepth 1 -maxdepth 1 -type d > $f_dirlist";

$datestring = gmtime();
print "Validation started at $datestring...\n";

my $cflag=0;
while($cflag!=1){
  my $getdirlist = "find $parentdir/output/glf_depth/chr$chr -mindepth 1 -maxdepth 1 -type d > $f_dirlist";
  &forkExecWait($getdirlist);
  open my $dirlist, '<', $f_dirlist or die "can't open $f_dirlist: $!";

  while(<$dirlist>){
    chomp;

    if(exists($filehash{$_}) && $filehash{$_}!=1){
      # print "Validating files in $_...";
      my @filecount=glob("$_/*.dp");
      my $numfiles = scalar @filecount;
      # chomp($numfiles);
      my $chknum=`wc -l $_/samples.ok | cut -d" " -f1`;
      chomp($chknum);

      if($chknum==$numind){
        my $perlcmd = "perl $parentdir/smaug-genetics/data_mgmt/process_glf_meandp.pl --chr $chr --dir $_";
        &forkExecWait($perlcmd);
        $filehash{$_}=1;

        my @path = split /\//, $_;
        my $range = $path[-1];
        my $meandp = "$parentdir/output/glf_depth/chr$chr.$range.txt";

        if(-e $meandp){
          print "Removing files in $_/\n";
          my $rmdpcmd="rm -f $_/\*.dp";
          &forkExecWait($rmdpcmd);
          my $rmokcmd="rm -f $_/samples.ok";
          &forkExecWait($rmokcmd);
        }
      }
    } else {
      $filehash{$_}=0;
    }
  }

  my $logfile="$parentdir/output/glf_depth/chr$chr.data.log";
  my $logcmd="sacct -j $ID --format=JobID,State | awk 'NR>2 {print \$2}' | sort | uniq | paste -d- -s > $logfile";
  &forkExecWait($logcmd);
  open my $log, '<', $logfile or die "can't open $logfile: $!";

  while(<$log>){
    chomp;
    if($_ eq "COMPLETED"){
      my @dirs;
      while(<$dirlist>){
        chomp;
        push @dirs, $_;
      }
      if(is_folder_empty($dirs[-1])){
        $cflag=1;
        print "VALIDATION COMPLETE\n";
        last;
      }
    }
  }

  sleep 120;
}

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
