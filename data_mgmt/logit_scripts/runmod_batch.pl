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

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait);

my $jobids = "1-4096";
my $mem=4000;

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
		$mem = 8000;
}

# foreach my $categ (@categs){
  # my $jobcmd="${categ}_logmod";
  my $builddatbatch = "$parentdir/slurm/logmod.txt";
  open my $mdFH, '>', $builddatbatch or die "can't write to $builddatbatch: $!\n";
  print $mdFH "#!/bin/bash \n";
	print $mdFH "#SGRIDBATCH CAT='AT' 'GC'\n";
	print $mdFH "#SGRIDBATCH INDEX=\$(seq 1 4096)\n";
  print $mdFH "#SBATCH --mail-type=FAIL \n";
  print $mdFH "#SBATCH --mail-user=$email \n";
  print $mdFH "#SBATCH --ntasks=1 \n";
  print $mdFH "#SBATCH --mem=8000 \n";
  print $mdFH "#SBATCH --time 00:20:00 \n";
  # print $mdFH "#SBATCH --job-name=logmod \n";
  print $mdFH "#SBATCH --partition=bipolar \n";
  # print $mdFH "#SBATCH --array=$jobids \n";
  print $mdFH "#SBATCH --requeue \n";
	# print $mdFH "#SBATCH --exclude=hunt-mc05,hunt-mc06,hunt-mc07,hunt-mc08,twins-mc01,twins-mc04,finnseq-mc02,dl3614,dl3615,dl3616,dl3617,dl3618,dl3619 \n";
  # print $mdFH "#SBATCH --exclude=psoriasis-mc01,psoriasis-mc02 \n";
  # print $mdFH "#SBATCH --output=\"$slurmdir/slurmJob-%J.out\" --error=\"$slurmdir/slurmJob-%J.err\" \n";
	print $mdFH "export STDOUT=\"$slurmdir/\$CAT.\$INDEX.out\"\n";
	print $mdFH "export STDERR=\"$slurmdir/\$CAT.\$INDEX.err\"\n";

	print $mdFH "if [ -f \$STDOUT ] ; then \n";
  print $mdFH "echo results for job \$STDOUT already exists \n";
  print $mdFH "exit 0 \n";
	print $mdFH "fi \n";

  print $mdFH "srun Rscript $parentdir/smaug-genetics/R/log_mod.r \$CAT $parentdir $libpath \$INDEX 1>\$STDOUT 2>\$STDERR \n";
  close($mdFH) or die "Unable to close file: $builddatbatch $!";

  my $slurmcmd="gridbatch $builddatbatch";
  forkExecWait($slurmcmd);
# }
