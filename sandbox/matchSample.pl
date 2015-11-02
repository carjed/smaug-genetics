#!/usr/local/bin/perl

##############################################################################
# Script gets all singletons which are adjacent to another singleton in the
# same individual
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

# Read in file with positions to match
my @POS;

my $f_positions = "/net/bipolar/jedidiah/chr20sing.singletons"; #main line for full processing
open my $positions, '<', $f_positions or die "can't open $f_positions: $!";

readline($positions);

my $outfile = "/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/adj_test.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

my $prevpos=0;
my $previous_line = "0\t0\tS\tA\tA";
my @prevA = split(/\t/, $previous_line);

my $testval=$prevA[4];

print "$testval\n";

while (<$positions>) {
	my @line=split(/\t/, $_);
	
	my $pos=$line[1];
	my $ind=$line[4];
	
	my @prevA = split(/\t/, $previous_line);
	
	if(($pos-1==$prevpos) & ($ind eq $prevA[4])){
	
		print OUT "$previous_line";
		print OUT "$_";
	
	}
	
	$previous_line=$_;
	$prevpos=$pos;
}



