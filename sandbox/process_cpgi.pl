#!/usr/local/bin/perl

##############################################################################
# Script to annotate singleton summary file with inside/outside CpGI status
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

my $chr;
GetOptions ('chr=i'=> \$chr);

## READ SINGLETON/COMMON DATA
my $chr10_sing = "/net/bipolar/jedidiah/mutation/output/chr$chr.expanded.summary";
open my $summ_sing, '<', $chr10_sing or die "can't open $chr10_sing: $!";
readline($summ_sing); #<-throws out summary header if it exists

# my $chr10_comm = "/net/bipolar/jedidiah/smaug-genetics/chr10.expanded.summary"; 
# open my $summ_comm, '<', $chr10_comm or die "can't open $chr10_comm: $!";
# readline($summ_comm); #<-throws out summary header if it exists

## DEFINE OUTPUT FILES
my $outfile_S = "chr$chr.cpgi.expanded.summary";
open(OUT_S, '>', $outfile_S) or die "can't write to $outfile_S: $!\n";

# my $outfile_C = "chr10.cpgi.expanded.summary";
# open(OUT_C, '>', $outfile_C) or die "can't write to $outfile_C: $!\n";

## GET CPG ISLAND DATA
my @cpgi_index;
my $f_cpgi = "/net/bipolar/jedidiah/reference_data/model-based-cpg-islands-hg19.txt";
open my $cpgi, '<', $f_cpgi or die "can't open $f_cpgi: $!";
my $header=readline($cpgi); #<-throws out summary header if it exists

# my $chr=10;

my $chrst;
while (<$cpgi>) {
	my @line=split(/\t/, $_);
		$chrst=$line[0];
		if ($chrst eq "chr$chr") {
			push(@cpgi_index, $_);
		}	
}

## INITIALIZE ARRAYS FOR ROWS AND POSITIONS
my @POS_S;
my @NEWSUMM_S;
my @POS_C;
my @NEWSUMM_C;
my $a_nu_start_cpg=0;


while (<$summ_sing>) {
	push (@POS_S, (split(/\t/, $_))[1]);
	push (@NEWSUMM_S, $_);
}

print OUT_S "$header\n";
foreach my $row (@NEWSUMM_S) {
	chomp $row;
	my @line=split(/\t/, $row);
	my $pos=$line[1];
	my $cpgi_ind = &getCpGI($pos);

	print OUT_S "$row\t$cpgi_ind\n";
}

# $a_nu_start_cpg=0;

# while (<$summ_comm>) {
	# push (@POS_C, (split(/\t/, $_))[1]);
	# push (@NEWSUMM_C, $_);
# }

# foreach my $row (@NEWSUMM_C) {
	# chomp $row;
	# my @line=split(/\t/, $row);
	# my $pos=$line[1];
	# my $cpgi = &getCpGI($pos);

	# print OUT_C "$row\t$cpgi\n";
# }


sub getCpGI {
	my $site = shift;
	my $hit=0;	
	
	my @first = split(/\t/, $cpgi_index[0]);
	# print "$first[1]\n";
	my @last = split(/\t/, $cpgi_index[$#cpgi_index]);	
	
	if ($site < $first[1]) {
		return $hit;
		last;
	} elsif ($site > $last[2]) {
		return $hit;
		last;	
	}	
	
	for my $i ($a_nu_start_cpg .. $#cpgi_index) {
		my $cpgi_int = $cpgi_index[$i];
		my @pair = split(/\t/, $cpgi_int);
		
		if (($site >= $pair[1]) && ($site <= $pair[2])) {
			$hit=1;
			# print "Site: $site is in CpG island: $pair[1], $pair[2]\n";
			$a_nu_start_cpg=$i;
			return $hit;
			last;
		}
	}
	
	return $hit;
}