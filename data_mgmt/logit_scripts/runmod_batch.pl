#!/usr/local/bin/perl

##############################################################################
# Script builds slurm batch file with specified parameters
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

my $email = $config->{email};
my $parentdir = $config->{parentdir};
my $libpath = $config->{libpath};

my $today = POSIX::strftime('%Y%m%d', localtime);
my $slurmdir = "$parentdir/output/slurm/$today";
  make_path("$slurmdir");

my $parentjob=1;
GetOptions('parentjob=i' => \$parentjob);

my $jobids = "1-4096";

# can run with --parentjob {jobno} to query failed child jobs and rebuild batch
# my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );
my @categs = qw(AT GC);
if($parentjob>1){
		$jobids=`sacct -j $parentjob --format=jobid%30,jobname%30,state | grep "FAILED" | grep -v "\\\." | awk '{print \$1}' | sed 's/${parentjob}_//g'`;
		$jobids =~ s/\r?\n/,/g;
		$jobids =~ s/,$//; # get rid of last comma

		# In re-run, only want loop below to run for category with failed jobs
		my $jobname=`sacct -j $parentjob --format=jobid%30,jobname%30,state | grep "FAILED" | grep -v "\\\." | awk 'NR==3{print \$2}'`;
		chomp $jobname;
		my $cat = substr($jobname, 0, 5);
		my @categs = ($cat);
}

foreach my $categ (@categs){
  my $jobcmd="${categ}_logmod";
  my $builddatbatch = "$parentdir/slurm/$jobcmd.txt";
  open my $mdFH, '>', $builddatbatch or die "can't write to $builddatbatch: $!\n";
  print $mdFH "#!/bin/sh \n";
  print $mdFH "#SBATCH --mail-type=FAIL \n";
  print $mdFH "#SBATCH --mail-user=$email \n";
  print $mdFH "#SBATCH --ntasks=1 \n";
  print $mdFH "#SBATCH --mem=6000 \n";
  print $mdFH "#SBATCH --time 10:00:00 \n";
  print $mdFH "#SBATCH --job-name=$jobcmd \n";
  print $mdFH "#SBATCH --partition=nomosix \n";
  print $mdFH "#SBATCH --array=$jobids \n";
  print $mdFH "#SBATCH --requeue \n";
	print $mdFH "#SBATCH --exclude=hunt-mc05,hunt-mc06,hunt-mc07,hunt-mc08,dl3614,dl3615,dl3616,dl3617,dl3618,dl3619 \n";
  # print $mdFH "#SBATCH --exclude=psoriasis-mc01,psoriasis-mc02 \n";
  print $mdFH "#SBATCH --output=\"$slurmdir/slurmJob-%J.out\" --error=\"$slurmdir/slurmJob-%J.err\" \n";
  print $mdFH "srun Rscript $parentdir/smaug-genetics/R/log_mod.r $categ $parentdir $libpath \$SLURM_ARRAY_TASK_ID\n";
  close($mdFH) or die "Unable to close file: $builddatbatch $!";

  my $slurmcmd="sbatch $builddatbatch";
  forkExecWait($slurmcmd);
}
