#!/usr/local/bin/perl

##############################################################################
# Script sums predicted rates over fixed-width windows
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

my $chr;
my $cat='';
my $bink='';

GetOptions('chr=s' => \$chr,
			'cat=s' => \$cat,
			'bink=s' => \$bink);

my $nextchr='';
if ($chr lt '22') {
	$nextchr=$chr+1;
} elsif ($chr eq '22') {
	$nextchr="X";
} else {
	$nextchr="Y";
}

my $f_data = "$parentdir/output/predicted/chr${chr}_${cat}_predicted.txt";
my $outfile = "$parentdir/output/predicted/chr${chr}_${cat}_binned.txt";
my $cpgoutfile = "$parentdir/output/predicted/chr${chr}_${cat}_binned_cpg.txt";

# initialize input file
open my $data, '<', $f_data or die "can't open $f_data: $!";

# initialize output
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";
open(OUTCPG, '>', $cpgoutfile) or die "can't write to $cpgoutfile: $!\n";

my $seq=&getRef();
my $altseq=$seq;
$altseq =~ tr/ACGT/TGCA/;

my $seqlength=length($seq);

# print "Hashing 1Mb bins...\n";
# my %hash=();
# while (<$coefs>){
# 	chomp;
# 	my @line=split(/\t/, $_);
# 	my $key=$line[0];
# 	my @betas=@line[1..$#line];
#
# 	$hash{$key}=[@betas];
# }

# <$data> for (1 .. 976_000_000);
my $adj=1;
# my $subseq=3;

my $firstLine = <$data>;
my @linearr=split(/\t/, $firstLine);
my $PREVBIN = ceil($linearr[1]/10);
my $SUM=0;
my $CPGSUM=0;
seek $data, 0, 0;
while (<$data>){
	chomp;
	my @line=split(/\t/, $_);
	my $POS=$line[0];
	my $BIN=ceil($line[1]/10); # Bin no. in predicted data is 100kb; coerce to 1Mb
	my $MU=$line[2];

	my $localseq = substr($seq, $POS-1, 2);
	my $altlocalseq = reverse substr($altseq, $POS-1, 2);

	if(($localseq eq "CG") | ($altlocalseq eq "CG")){
		if($BIN==$PREVBIN){
			$CPGSUM+=$MU;
		} else {
			print OUTCPG "$chr\t$PREVBIN\t$CPGSUM\n";
			$CPGSUM=0;
			$PREVBIN = $BIN;
		}
	} else {
		if($BIN==$PREVBIN){
			$SUM+=$MU;
		} else {
			print OUT "$chr\t$PREVBIN\t$SUM\n";
			$SUM=0;
			$PREVBIN = $BIN;
		}
	}
}

sub getRef{
	my $f_fasta="$parentdir/reference_data/human_g1k_v37.fasta";
	my $mask_flag=0;
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
	if($mask_flag==1){
		while (<$fasta>) {
			chomp;
			if (/^>$chr$/../^>$nextchr$/) {
				next if /^>$chr$/ || /^>$nextchr$/;
				$seq .=$_;
			}
		}
	} else {
		while (<$fasta>) {
			chomp;
			if (/>$chr /../>$nextchr /) {
				next if />$chr / || />$nextchr /;
				$seq .=$_;
			}
		}
	}

	return $seq;
}
