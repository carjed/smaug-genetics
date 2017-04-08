#!/usr/local/bin/perl

##############################################################################
# Sort predicted rates
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

my $adj = 3;
my $data = $config->{data};
my $email = $config->{email};
my $parentdir = $config->{parentdir};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

my $today = POSIX::strftime('%Y%m%d', localtime);
my $slurmdir = "$parentdir/output/slurm/$today";
  make_path("$slurmdir");

my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );

foreach my $categ (@categs){
  my $jobcmd="${categ}_sort_pred";
  my $builddatbatch = "$parentdir/slurm/$jobcmd.txt";
  open my $mdFH, '>', $builddatbatch or die "can't write to $builddatbatch: $!\n";
  print $mdFH "#!/bin/sh \n";
  print $mdFH "#SBATCH --mail-type=FAIL \n";
  print $mdFH "#SBATCH --mail-user=$email \n";
  print $mdFH "#SBATCH --ntasks=1 \n";
  print $mdFH "#SBATCH --mem=4000 \n";
  print $mdFH "#SBATCH --time 2:00:00 \n";
  print $mdFH "#SBATCH --job-name=$jobcmd \n";
  print $mdFH "#SBATCH --partition=nomosix \n";
  print $mdFH "#SBATCH --array=1-22 \n";
  print $mdFH "#SBATCH --requeue \n";
  # print $mdFH "#SBATCH --exclude=psoriasis-mc01,psoriasis-mc02 \n";
  print $mdFH "#SBATCH --output=\"$slurmdir/slurmJob-%J.out\" --error=\"$slurmdir/slurmJob-%J.err\" \n";
  print $mdFH "srun cat $parentdir/output/predicted/$categ/chr\${SLURM_ARRAY_TASK_ID}/* | sort -V -k2 > $parentdir/output/predicted/chr\${SLURM_ARRAY_TASK_ID}.$categ.txt";
  close($mdFH) or die "Unable to close file: $builddatbatch $!";

  my $slurmcmd="sbatch $builddatbatch";
  forkExecWait($slurmcmd);
}
