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

# foreach my $categ (@categs){
  # my $jobcmd="${categ}_sort_pred";
  my $builddatbatch = "$parentdir/slurm/sort_pred_batch.txt";
  open my $mdFH, '>', $builddatbatch or die "can't write to $builddatbatch: $!\n";
  print $mdFH "#!/bin/bash \n";
  print $mdFH "#SGRIDBATCH CAT='AT_CG' 'AT_GC' 'AT_TA' 'GC_AT' 'GC_CG' 'GC_TA'\n";
  print $mdFH "#SGRIDBATCH INDEX=\$(seq 1 22)\n";
  print $mdFH "#SBATCH --mail-type=FAIL \n";
  print $mdFH "#SBATCH --mail-user=$email \n";
  print $mdFH "#SBATCH --ntasks=1 \n";
  print $mdFH "#SBATCH --mem=4000 \n";
  print $mdFH "#SBATCH --time 2:00:00 \n";
  print $mdFH "#SBATCH --job-name=sort_pred \n";
  print $mdFH "#SBATCH --partition=bipolar \n";
  # print $mdFH "#SBATCH --array=1-22 \n";
  print $mdFH "#SBATCH --requeue \n";
  print $mdFH "export STDOUT=\"$slurmdir/chr\$INDEX.\$CAT.out\"\n";
  print $mdFH "export STDERR=\"$slurmdir/chr\$INDEX.\$CAT.err\"\n";
  # print $mdFH "#SBATCH --exclude=psoriasis-mc01,psoriasis-mc02 \n";
  # print $mdFH "#SBATCH --output=\"$slurmdir/slurmJob-%J.out\" --error=\"$slurmdir/slurmJob-%J.err\" \n";
  print $mdFH "srun cat $parentdir/output/predicted/\$CAT/chr\$INDEX/* | sort -V -k2 > $parentdir/output/predicted/chr\$INDEX.\$CAT.txt";
  close($mdFH) or die "Unable to close file: $builddatbatch $!";

  my $slurmcmd="gridbatch $builddatbatch";
  forkExecWait($slurmcmd);
# }
