#!/usr/local/bin/perl

##############################################################################
# Creates index file for GotCloud to run variant calling on trio bams
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
foreach (@samples) {
	print OUT "$_\t$cov\t/net/bipolar/jedidiah/sardinia/${cov}x-in/${_}_${cov}x.bam\n";
}