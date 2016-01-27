#!/usr/local/bin/perl

##############################################################################
# Downsampling pipeline for trio data
#
# Runs the following sequence of commands:
# 1. Downsample full BAM files to desired coverage
# 2. Index resulting BAMs
# 3. Write gotCloud index file
# 4. Run gotCloud snpcall and ldrefine
# 5. Merge all output chromosome vcfs to single file
# 6. Count per-subject singletons
#
# Currently must run for each downsampled coverage level
##############################################################################

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Cwd;
use Benchmark;

my $cov;

GetOptions ('cov=i'=> \$cov);

my $prop=$cov/80;
my @samples = ("28009", "28038", "28047");

#print "$prop\n";

my $outfile = "/net/bipolar/jedidiah/gotcloudExample/GBR60bam.index";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

#subsample BAMs and index; output index file for gotCloud
my $downsampleCMD;
foreach (@samples) {
	print OUT "SD$_\t$cov\t/net/bipolar/jedidiah/sardinia/${cov}x-in/${_}_${cov}x.bam\n";

	$downsampleCMD = "/net/wonderland/home/yancylo/bin/samtools-0.2.0-rc10/samtools view -uh -o /net/bipolar/jedidiah/sardinia/${cov}x-in/${_}_${cov}x.bam -s $prop /net/sardinia/progenia/SardiNIA/paraparesi/80x/$_.merged80x.bam";
	&forkExecWait($downsampleCMD);
}

my $indexCMD;	
foreach (@samples) {
	$indexCMD = "/usr/cluster/bin/samtools index /net/bipolar/jedidiah/sardinia/${cov}x-in/${_}_${cov}x.bam";
	&forkExecWait($indexCMD);
	
}

# foreach (@samples) {
	# print OUT "$_\t$cov\t/net/bipolar/jedidiah/sardinia/${cov}x-in/${_}_${cov}x.bam\n";
# }

#gotCloud snpcall and ldrefine
my $snpcallCMD = "/usr/cluster/bin/gotcloud snpcall --conf /net/bipolar/jedidiah/gotcloudExample/GBR60vc.conf --outdir /net/bipolar/jedidiah/sardinia/${cov}x-gw --numjobs 10";
&forkExecWait($snpcallCMD);

my $ldrefineCMD = "/usr/cluster/bin/gotcloud ldrefine --conf /net/bipolar/jedidiah/gotcloudExample/GBR60vc.conf --outdir /net/bipolar/jedidiah/sardinia/${cov}x-gw --numjobs 10";
&forkExecWait($ldrefineCMD);

#merge all chromosomes to single VCF
my @vcfs;
my $root="/net/bipolar/jedidiah/sardinia/${cov}x-gw/beagle";

for my $i (1 .. 22) {
	push (@vcfs, "$root/chr$i/chr$i.filtered.PASS.beagled.vcf.gz");
}

my $input;
foreach (@vcfs) {
	#chomp;
	$input .= "$_ ";
}

my $concatCMD = "/net/bipolar/jedidiah/bcftools/bcftools concat -O z $input > /net/bipolar/jedidiah/sardinia/${cov}x_full.vcf.gz";
&forkExecWait($concatCMD);

#Count per-subject singletons
my $singletonCMD= "/usr/cluster/bin/vcftools --gzvcf /net/bipolar/jedidiah/sardinia/${cov}x_full.vcf.gz --singletons --stdout | cut -f5 | sort | uniq -c >> ${cov}x.singletons";
&forkExecWait($singletonCMD);

##############################################################################
# fork-exec-wait subroutine
##############################################################################
sub forkExecWait {
    my $cmd = shift;
    #print "forkExecWait(): $cmd\n";
    my $kidpid;
    if ( !defined($kidpid = fork()) ) {
		die "Cannot fork: $!";
    }
    elsif ( $kidpid == 0 ) {
		exec($cmd);
		die "Cannot exec $cmd: $!";
    }
    else {
		waitpid($kidpid,0);
    }
}