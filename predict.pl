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
my $binw='';

GetOptions('chr' => \$chr,
			'cat' => \$cat,
			'binw' => \$binw);
			
my $f_coefs = "$parentdir/output/logmod_data/${cat}_${binw}_coefs.txt";
my $f_data = "$parentdir/output/logmod_data/chr${chr}_${cat}_sites.txt";
my $outfile = "$parentdir/output/predicted/chr${chr}_${cat}_predicted.txt";

# initialize covariate data
open my $coefs, '<', $f_coefs or die "can't open $f_coefs: $!";

# initialize input file
open my $data, '<', $f_data or die "can't open $f_data: $!";

# initialize output
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";


print "Hashing model parameters...\n";
my %hash=();
while (<$coefs>){
	chomp;
	my @line=split(/\t/, $_);
	my $key=$line[0];
	my @betas=@line[1..$#line];
	
	$hash{$key}=[@betas];
}

print "Getting predicted values...\n";
my $rownum=0;
my $mbnum=0;
my $sum=0;
my $mean;

# <$data> for (1 .. 976_000_000);
while (<$data>){
	chomp;
	my @line=split(/\t/, $_);
	my $CHR=$line[0];
	my $POS=$line[2];
	my $SEQ=$line[3];
	
	if($SEQ !~ /[MNSW]/){
		my @vals=(1,@line[5..$#line]);
		my @betas=@{$hash{$SEQ}};
		
		# print "$vals[$_]\n" for 0 .. $#vals;
		# print "$betas[$_]\n" for 0 .. $#betas;
		
		my $pred=0;
		$pred+=$vals[$_]*$betas[$_] for 0 .. $#vals;
		
		$pred=nearest(0.001, ((exp($pred)/(exp($pred)+1))/314244)*1e8);
		# $pred=nearest(0.001, ((exp($pred)/(exp($pred)+1))/314244)*1e8);
		
		print OUT "$CHR\t$POS\t$pred\n";
		
		$rownum++;
		$sum+=$pred;
		$mean=$sum/$rownum;
		
		if($rownum%10000000==0){
			$mbnum+=10;
			print "Finished $mbnum million bases (current mean rate: $mean)...\n";
		}
	}
}