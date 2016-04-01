#!/usr/local/bin/perl

##############################################################################
# Script reads full data file and matrix of model coefficients for given
# mutation category and outputs per-site predicted mutation rates
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

my $f_data = "$parentdir/output/predicted/chr${chr}_${cat}_predicted.txt";
my $outfile = "$parentdir/output/predicted/chr${chr}_${cat}_binned.txt";

# initialize input file
open my $data, '<', $f_data or die "can't open $f_data: $!";

# initialize output
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";
# 
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
my $firstLine = <$data>;
my @linearr=split(/\t/, $firstLine);
my $PREVBIN = $linearr[1];
my $SUM=0;
seek $data, 0, 0;
while (<$data>){
	chomp;
	my @line=split(/\t/, $_);
	my $POS=$line[0];
	my $BIN=$line[1];
	my $MU=$line[2];

	if($BIN==$PREVBIN){
		$SUM+=$MU;
	} else {
		print OUT "$chr\t$BIN\t$SUM\n";
		$SUM=0;
	}

	my $PREVBIN = $BIN;

}
