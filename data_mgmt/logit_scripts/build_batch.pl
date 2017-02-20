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

my $trainchr='';
my $cat='';
my $bink='';

GetOptions('trchr=s' => \$trainchr,
			'cat=s' => \$cat,
			'bink=s' => \$bink);

print "$trainchr\n";
print "$cat\n";

# initialize output
my $outfile = "$parentdir/smaug-genetics/slurm_predict.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

print OUT "#!/bin/sh \n";
print OUT "#SBATCH --mail-type=FAIL \n";
print OUT "#SBATCH --mail-user=jedidiah\@umich.edu \n";
print OUT "#SBATCH --ntasks=1 \n";
print OUT "#SBATCH --mem=6000 \n";
print OUT "#SBATCH --time 08:00:00 \n";
print OUT "#SBATCH --job-name=predict_$cat \n";
print OUT "#SBATCH --partition=nomosix \n";
print OUT "#SBATCH --array=$trainchr \n";
print OUT "#SBATCH --exclude=topmed,topmed2 \n";
print OUT "#SBATCH --output=\"/net/bipolar/jedidiah/mutation/output/slurm/slurmJob-%J.out\" --error=\"/net/bipolar/jedidiah/mutation/output/slurm/slurmJob-%J.err\" \n";
print OUT "srun perl /net/bipolar/jedidiah/mutation/smaug-genetics/predict.pl --chr \${SLURM_ARRAY_TASK_ID} --cat $cat --bink $bink \n";
