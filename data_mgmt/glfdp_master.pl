#!/usr/local/bin/perl

##############################################################################
# Scans summarized glf files (every 10bp; colnames: chr, pos, ref, dp) and
# outputs mean depth
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
use Cwd;
use Benchmark;
use Tie::File;

my $chr;
my $help=0;
my $man=0;

# Set options and inputs
GetOptions ('chr=i'=> \$chr,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;
# my $chr=1;

my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";

my %filehash=();
my $f_dirlist = "$parentdir/output/glf_depth/chr${chr}_glf_dirlist.txt";
my $getdirlist = "find $parentdir/output/glf_depth/chr$chr -mindepth 1 -maxdepth 1 -type d > $f_dirlist";
&forkExecWait($getdirlist);
open my $dirlist, '<', $f_dirlist or die "can't open $f_dirlist: $!";

while(<$dirlist>){
  chomp;

  if(exists($filehash{$_}) && $filehash{$_}!=1){
    # print "Validating files in $_...";
    my @filecount=glob("$_/*.dp");
    my $numfiles = scalar @filecount;
    # chomp($numfiles);
    my $chknum=`wc -l $_/samples.ok | cut -d" " -f1`;
    chomp($chknum);

    if($chknum==$numind){
      my $perlcmd = "perl $parentdir/smaug-genetics/data_mgmt/process_glf_meandp.pl --chr $chr --dir $_";
      &forkExecWait($perlcmd);
      $filehash{$_}=1;

      my @path = split /\//, $_;
      my $range = $path[-1];
      my $meandp = "$parentdir/output/glf_depth/chr$chr.$range.txt";

      if(-e $meandp){
        print "Removing files in $_/\n";
        my $rmdpcmd="rm -f $_/\*.dp";
        &forkExecWait($rmdpcmd);
        my $rmokcmd="rm -f $_/samples.ok";
        &forkExecWait($rmokcmd);
      }
    }
  } else {
    $filehash{$_}=0;
  }
}

# initialize and run sbatch file
my $outfile = "$parentdir/smaug-genetics/data_mgmt/slurm_glf_meandp.$chr.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";
print OUT "#!/bin/sh \n";
print OUT "#SBATCH --mail-type=FAIL \n";
print OUT "#SBATCH --mail-user=jedidiah\@umich.edu \n";
print OUT "#SBATCH --ntasks=1 \n";
print OUT "#SBATCH --mem=2000 \n";
print OUT "#SBATCH --time 2:00:00 \n";
print OUT "#SBATCH --job-name=chr${chr}_glf_meandp \n";
print OUT "#SBATCH --partition=nomosix \n";
print OUT "#SBATCH --array=1-$numjobs \n"; # change to 1-$numjobs
print OUT "#SBATCH --output=\"$parentdir/output/slurm/slurmJob-%J.out\" --error=\"$parentdir/output/slurm/slurmJob-%J.err\" \n";
print OUT "srun perl $parentdir/smaug-genetics/data_mgmt/process_glf_meandp.pl --chr $chr --ind \${SLURM_ARRAY_TASK_ID}";
close(OUT) or die "Unable to close file: $outfile $!";

my $slurmcmd="sbatch $outfile";
&forkExecWait($slurmcmd);
