#!/usr/local/bin/perl

##############################################################################
# Script builds slurm batch file to add depth column for logit model
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
my $today = POSIX::strftime('%Y%m%d', localtime);
my $slurmdir = "$parentdir/output/slurm/$today";
  make_path("$slurmdir");

my $numjobs=4096;
my $categ;

my $help=0;
my $man=0;

# Set options and inputs
GetOptions ('categ=s'=> \$categ,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

my $jobcmd="${categ}_add_dp";
my $workerbatch = "$parentdir/slurm/slurm_$jobcmd.txt";
open my $wFH, '>', $workerbatch or die "can't write to $workerbatch: $!\n";
print $wFH "#!/bin/sh \n";
print $wFH "#SBATCH --mail-type=FAIL \n";
print $wFH "#SBATCH --mail-user=jedidiah\@umich.edu \n";
print $wFH "#SBATCH --ntasks=1 \n";
print $wFH "#SBATCH --mem=6000 \n";
print $wFH "#SBATCH --time 20:00:00 \n";
print $wFH "#SBATCH --job-name=$jobcmd \n";
print $wFH "#SBATCH --partition=nomosix \n";
print $wFH "#SBATCH --array=1-$numjobs \n"; # change to 1-$numjobs
print $wFH "#SBATCH --requeue \n";
print $wFH "#SBATCH --exclude=inpsyght \n";
# print $wFH "#SBATCH --exclude=psoriasis-mc01,psoriasis-mc02 \n";
print $wFH "#SBATCH --output=\"$slurmdir/slurmJob-%J.out\" --error=\"$slurmdir/slurmJob-%J.err\" \n";
print $wFH "srun perl $parentdir/smaug-genetics/data_mgmt/add_dp_worker.pl --in \${SLURM_ARRAY_TASK_ID} --categ $categ \n";
close($wFH) or die "Unable to close file: $workerbatch $!";
my $slurmcmd="sbatch $workerbatch";
&forkExecWait($slurmcmd);

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
