#!/usr/local/bin/perl

##############################################################################
# Script used by obsexp6.r to obtain negative (i.e., non-mutated) examples
# used in logistic regression
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

# Set options and inputs
my $wdir=getcwd;
my $parentdir=dirname($wdir);

my $mask_flag=1;
my $adj=2;
my $binwidth=100000;

my $chr=20;

my $nextchr;
if ($chr<22) {
	$nextchr=$chr+1;
} elsif ($chr==22) {
	$nextchr="X";
}

my $subseq=2;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

my $outfile = "/net/bipolar/jedidiah/mutation/smaug-genetics/negative_examples.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

# my $seq="ACGTAAGCCGTAAC";
# 2,3,7,8,9,10,14

my $seq=&getRef();
my $altseq=$seq;
$altseq =~ tr/ACGT/TGCA/;

my $min=5;
my $max=length($seq)-5;

print "seqlength: $max\n";

print OUT "CHR \t POS \t BIN \t Sequence \n";

# Read in file with positions to exclude from analyses
my @POS;

my $f_positions = "/net/bipolar/jedidiah/mutation/smaug-genetics/exclusion_list.txt"; #main line for full processing
open my $positions, '<', $f_positions or die "can't open $f_positions: $!";

while (<$positions>) {
	next unless $_ =~ /^[^,]*$/;
	push (@POS, $_);
}

for my $i (1 .. 1000) {

	my $pos = $min + int(rand($max - $min));
	my $base = substr($seq, $pos-1, 1);

	if(($base =~ /G|C/) & !(grep(/^$pos$/, @POS))){
		push (@POS, $pos); # add position to exclusion list
		my $localseq = substr($seq, $pos-$adj-1, $subseq);
		my $altlocalseq = reverse substr($altseq, $pos-$adj-1, $subseq);
		my $bin = ceil($pos/$binwidth);
		
		# Coerce local sequence info to format used in R
		my $sequence;
		if(substr($localseq,$adj,1) lt substr($altlocalseq,$adj,1)){
			$sequence = $localseq . '(' . $altlocalseq . ')';
		} else {
			$sequence = $altlocalseq . '(' . $localseq . ')';
		}
		if ($sequence !~ /N/) {
			print OUT "$chr\t$pos\t$bin\t$sequence\n";
		} 
	}
}


sub getRef{
	my $f_fasta;
	if($mask_flag==1){
		$f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
	} else {
		$f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
	}

	if (-e $f_fasta) {
		print "Using reference genome: $f_fasta\n";
	} else {
		print "Reference genome not found in parent directory. Would you like to download one? (y/n): ";
		my $choice = <>;
		chomp $choice;
		if ($choice eq "y") {
			my $dlcmd="wget -P $parentdir/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz";
			&forkExecWait($dlcmd);
			my $unzipcmd="gunzip $parentdir/human_g1k_v37.fasta";
			&forkExecWait($unzipcmd);
		} else {
			die "Please upload an appropriate reference genome to the parent directory\n";
		}
	}

	open my $fasta, '<', $f_fasta or die "can't open $f_fasta: $!";

	##############################################################################
	# Retrieve reference sequence for selected chromosome
	# -also returns symmetric sequence to be used in local sequence analysis
	##############################################################################

	print "Getting reference sequence for chromosome $chr...\n";

	my $seq;
	while (<$fasta>) {
		chomp;
		if (/^>$chr$/../^>$nextchr$/) {
			next if /^>$chr$/ || /^>$nextchr$/;
			$seq .=$_;
		}
	}
	
	return $seq;
}