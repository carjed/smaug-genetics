#!/usr/local/bin/perl

##############################################################################
# Script sums predicted rates over fixed-width windows
# NEED TO FIX: removed bin column on 4/15/16
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
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));

my $config = LoadFile("$configpath/_config.yaml");

my $adj = $config->{adj};
my $mac = $config->{mac};
my $binw = $config->{binw};
my $data = $config->{data};
my $bin_scheme = $config->{bin_scheme};
my $parentdir = $config->{parentdir};
my $count_motifs = $config->{count_motifs};
my $expand_summ = $config->{expand_summ};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

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

my $f_data = "$parentdir/output/predicted/chr${chr}.${cat}.txt";
my $outfile = "$parentdir/output/predicted/binned/chrs/chr${chr}_${cat}_binned.txt";
# my $cpgoutfile = "$parentdir/output/predicted/chr${chr}_${cat}_binned_cpg.txt";

# initialize input file
open my $data, '<', $f_data or die "can't open $f_data: $!";

# initialize output
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";
# open(OUTCPG, '>', $cpgoutfile) or die "can't write to $cpgoutfile: $!\n";
my $f_fasta;
if($data eq "mask"){
  $f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
} else {
  $f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
}

my $seq=getRef($f_fasta, $chr);
# my $altseq=$seq;
# $altseq =~ tr/ACGT/TGCA/;

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
my $PREVBIN = ceil($linearr[2]/10);
my $SUM=0;
my $CPGSUM=0;
seek $data, 0, 0;

my $i=1;
while (<$data>){
	# if($i>10){last;}
	chomp;
	my @line=split(/\t/, $_);
	my $POS=$line[1];
	my $BIN=ceil($line[2]/10); # Bin no. in predicted data is 100kb; coerce to 1Mb
	my $MU=$line[3];

	my $localseq = substr($seq, $POS-2, 3);
	my $altlocalseq = reverse $localseq;
	$altlocalseq =~ tr/ACGT/TGCA/;

		if($BIN==$PREVBIN){
			if((substr($localseq, 1, 2) eq "CG") || (substr($altlocalseq, 1, 2) eq "CG")){
				# print "$POS\t$localseq\t$altlocalseq\tCpG\n";
				$CPGSUM+=$MU;
			} else {
				$SUM+=$MU;
				# print "$POS\t$localseq\t$altlocalseq\n";
			}
		} else {
			print OUT "$chr\t$PREVBIN\t$SUM\t$CPGSUM\n";
			$SUM=0;
			$CPGSUM=0;
			$PREVBIN = $BIN;
		}
		$i++;
}
