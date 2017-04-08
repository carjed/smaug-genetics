#!/usr/local/bin/perl

##############################################################################
# Script sums predicted rates over fixed-width windows
##############################################################################
use strict;
use warnings;
use POSIX;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));

my $config = LoadFile("$configpath/_config.yaml");

my $binw = $config->{binw};
my $data = $config->{data};
my $bin_scheme = $config->{bin_scheme};
my $parentdir = $config->{parentdir};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

my $chr;
my $cat='';

GetOptions('chr=s' => \$chr,
			'cat=s' => \$cat);

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
my $seqlength=length($seq);


my $firstLine = <$data>;
my @linearr=split(/\t/, $firstLine);
my $PREVBIN = ceil($linearr[1]/$binw);
my $SUM=0;
my $CPGSUM=0;
seek $data, 0, 0;

my $i=1;
while (<$data>){
	# if($i>10){last;}
	chomp;
	my @line=split(/\t/, $_);
	my $POS=$line[1];
	my $MU=$line[2];
	my $BIN=ceil($line[1]/$binw);

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
