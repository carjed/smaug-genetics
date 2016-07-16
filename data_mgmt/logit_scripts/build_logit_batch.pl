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

my $parentdir="/net/bipolar/jedidiah/mutation";

my $jobids=`sacct -j 27532888 --format=jobid%30,jobname%30,state | grep "FAILED" | grep "test_logit" | awk 'NR>1 {print \$1}' | sed 's/27532888_//g'`;

# print "$jobids\n";
$jobids =~ s/\r?\n/,/g;
$jobids =~ s/,$//; # get rid of last comma

# initialize output
my $outfile = "$parentdir/smaug-genetics/R/logit_batch_redo.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

print OUT "#!/bin/sh \n";
print OUT "#SBATCH --mail-type=FAIL \n";
print OUT "#SBATCH --mail-user=jedidiah\@umich.edu \n";
print OUT "#SBATCH --ntasks=1 \n";
print OUT "#SBATCH --mem=6000 \n";
print OUT "#SBATCH --time 10:00:00 \n";
print OUT "#SBATCH --job-name=logit_slurm \n";
print OUT "#SBATCH --partition=nomosix \n";
print OUT "#SBATCH --array=$jobids \n";
print OUT "#SBATCH --requeue \n";
# print OUT "#SBATCH --exclude=topmed,topmed2 \n";
print OUT "#SBATCH --output=\"/net/bipolar/jedidiah/mutation/output/slurm/slurmJob-%J.out\" --error=\"/net/bipolar/jedidiah/mutation/output/slurm/slurmJob-%J.err\" \n";
print OUT "srun Rscript log_mod.r --categ=AT_GC --jobid=\$SLURM_ARRAY_TASK_ID\n";
