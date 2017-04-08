#!/usr/local/bin/perl

##############################################################################
# Used to obtain full data from logistic regression model
# loops through reference genome and outputs 1 line per base, as long as
# valid covariate data exists
##############################################################################

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils 'pairwise';
use Cwd;
use Benchmark;
use Tie::File;
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

my $adj = 3;
my $data = $config->{data};
my $parentdir = $config->{parentdir};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );

foreach my $categ (@categs){
  my $builddatcmd = "perl $relpath/getNonMut.pl --categ $categ --chr $chr ";
  forkExecWait($builddatcmd);

  my $builddatbatch = "$parentdir/slurm/builddat_$categ.txt";
  open my $mdFH, '>', $builddatbatch or die "can't write to $builddatbatch: $!\n";
  print $mdFH "#!/bin/sh \n";
  print $mdFH "#SBATCH --mail-type=FAIL \n";
  print $mdFH "#SBATCH --mail-user=$email \n";
  print $mdFH "#SBATCH --ntasks=1 \n";
  print $mdFH "#SBATCH --mem=2000 \n";
  print $mdFH "#SBATCH --time 20:00:00 \n";
  print $mdFH "#SBATCH --job-name=$jobcmd \n";
  print $mdFH "#SBATCH --partition=nomosix \n";
  print $mdFH "#SBATCH --array=1-22 \n";
  print $mdFH "#SBATCH --requeue \n";
  # print $mdFH "#SBATCH --exclude=psoriasis-mc01,psoriasis-mc02 \n";
  print $mdFH "#SBATCH --output=\"$slurmdir/slurmJob-%J.out\" --error=\"$slurmdir/slurmJob-%J.err\" \n";
  print $mdFH "srun perl $relpath/getNonMut.pl --categ $categ --chr \${SLURM_ARRAY_TASK_ID}";
  close($mdFH) or die "Unable to close file: $builddatbatch $!";

  $slurmcmd="sbatch $meandpbatch";
  forkExecWait($slurmcmd);
}


# my $fullfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_full.txt.gz";
#
# my $subcmd = "zcat $fullfile | awk '{print >> \"$parentdir/output/logmod_data/motifs/$categ/${categ}_\" \$4 \".txt\"}'";
# forkExecWait($subcmd);
