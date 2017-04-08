#!/usr/local/bin/perl

##############################################################################
# Concatenate per-motif predicted rates by chromosome and sort by position
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

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

my $today = POSIX::strftime('%Y%m%d', localtime);
my $slurmdir = "$parentdir/output/slurm/$today";
  make_path("$slurmdir");

my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );

foreach my $categ (@categs){
  my $jobcmd="${categ}_bin_pred";
  my $builddatbatch = "$parentdir/slurm/$jobcmd.txt";
  open my $mdFH, '>', $builddatbatch or die "can't write to $builddatbatch: $!\n";
  print $mdFH "#!/bin/sh \n";
  print $mdFH "#SBATCH --mail-type=FAIL \n";
  print $mdFH "#SBATCH --mail-user=$email \n";
  print $mdFH "#SBATCH --ntasks=1 \n";
  print $mdFH "#SBATCH --mem=1000 \n";
  print $mdFH "#SBATCH --time 00:30:00 \n";
  print $mdFH "#SBATCH --job-name=$jobcmd \n";
  print $mdFH "#SBATCH --partition=nomosix \n";
  print $mdFH "#SBATCH --array=1-22 \n";
  print $mdFH "#SBATCH --requeue \n";
  # print $mdFH "#SBATCH --exclude=psoriasis-mc01,psoriasis-mc02 \n";
  print $mdFH "#SBATCH --output=\"$slurmdir/slurmJob-%J.out\" --error=\"$slurmdir/slurmJob-%J.err\" \n";
  print $mdFH "srun perl $relpath/bin_pred_worker.pl --chr \${SLURM_ARRAY_TASK_ID} --cat $categ";
  close($mdFH) or die "Unable to close file: $builddatbatch $!";

  my $slurmcmd="sbatch $builddatbatch";
  forkExecWait($slurmcmd);
}
