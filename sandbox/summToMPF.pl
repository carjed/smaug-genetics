#!/usr/local/bin/perl

##############################################################################
# Script to convert genome-wide summary file to Mutation Position Format, 
# as used in pmsignature R package
##############################################################################

use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Path qw(make_path);

# Set options and inputs
my $wdir=getcwd;
my $parentdir=dirname($wdir);

my $outfile = "/net/bipolar/jedidiah/mutation/output/full_summ.mpf";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

my $summfile = "/net/bipolar/jedidiah/mutation/output/5bp_100k/full.summary"; #main line for full processing
open my $summ, '<', $summfile or die "can't open $summfile: $!";

readline($summ); #<-throws out summary header if it exists

while (<$summ>) {
	next unless $_ =~ /^[^,]*$/;
	my $cutoff = 4;
	my $newline = join("\t", (split(/\t/, $_, $cutoff+1))[0..$cutoff-1]);
	print OUT "1\t$newline\n";
}