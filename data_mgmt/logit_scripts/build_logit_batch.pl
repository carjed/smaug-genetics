#!/usr/local/bin/perl

##############################################################################
# Script builds slurm batch file with specified parameters
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
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

my $email = $config->{email};
my $parentdir = $config->{parentdir};
my $libpath = $config->{libpath};

my $categ="AT_GC";
my $parentjob=1;

GetOptions('categ=s' => \$categ,
			'parentjob=i' => \$parentjob);

my $jobids = "1-4096";

# can run with --parentjob {jobno} to query failed child jobs and rebuild batch
if($parentjob>1){
		$jobids=`sacct -j $jobid --format=jobid%30,jobname%30,state | grep "FAILED" | grep -v "\\\." | awk '{print \$1}' | sed 's/${parentjob}_//g'`;
		$jobids =~ s/\r?\n/,/g;
		$jobids =~ s/,$//; # get rid of last comma
}

# initialize output
my $outfile = "$parentdir/slurm/logit_batch_redo.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

print OUT "#!/bin/sh \n";
print OUT "#SBATCH --mail-type=FAIL \n";
print OUT "#SBATCH --mail-user=$email \n";
print OUT "#SBATCH --ntasks=1 \n";
print OUT "#SBATCH --mem=6000 \n";
print OUT "#SBATCH --time 10:00:00 \n";
print OUT "#SBATCH --job-name=logit_slurm_${categ} \n";
print OUT "#SBATCH --partition=nomosix \n";
print OUT "#SBATCH --array=$jobids \n";
print OUT "#SBATCH --requeue \n";
print OUT "#SBATCH --exclude=hunt-mc05,hunt-mc06,hunt-mc07,hunt-mc08,dl3614,dl3615,dl3616,dl3617,dl3618,dl3619 \n";
# print OUT "#SBATCH --exclude=topmed,topmed2 \n";
print OUT "#SBATCH --output=\"/net/bipolar/jedidiah/mutation/output/slurm/slurmJob-%J.out\" --error=\"/net/bipolar/jedidiah/mutation/output/slurm/slurmJob-%J.err\" \n";
print OUT "srun Rscript $parentdir/smaug-genetics/R/log_mod.r $categ $parentdir $libpath \$SLURM_ARRAY_TASK_ID\n";
