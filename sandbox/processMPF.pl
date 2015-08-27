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
use Tie::File;

# Set options and inputs
my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";

my $baseopt;
my $chr;
my $categ;
my $bw;

GetOptions ('b=s'=> \$baseopt,
'chr=s'=> \$chr,
'categ=s' => \$categ,
'bw=i' => \$bw) or pod2usage(1);

my $b1;
my $b2;
if($baseopt eq "AT"){
	$b1="A";
	$b2="T";
} else {
	$b1="C";
	$b2="G";
}

my $mask_flag=0;
my $adj=2;
my $binwidth=100000;

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

# initialize output file
my $outfile = "outdat.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

# initialize covariate data
my $f_covs = "$parentdir/output/logmod_data/${bw}kb_mut_cov2.txt";
open my $covs, '<', $f_covs or die "can't open $f_covs: $!";

# initialize input file
my $f_positions = "$parentdir/output/logmod_data/test.mpf";
open my $positions, '<', $f_positions or die "can't open $f_positions: $!";

# initialize phastCons data
my $f_cons = "$parentdir/reference_data/chr$chr.phastCons46way.primates.wigFix";
open my $cons, '<', $f_cons or die "can't open $f_cons: $!";

# Get reference sequence
my $seq=&getRef();
my $altseq=$seq;
$altseq =~ tr/ACGT/TGCA/;

# my $seqlength=length($seq);
# print "seqlength of chr$chr: $max\n"; #<-used to validate that getRef() returns correct seq length

my $printheader=0;
if($printheader==1){
	print OUT "CHR \t POS \t BIN \t Sequence \t mut \n"; #<-add header to output, if needed
}

# Create hash keyed by Chr/Bin pairs, with row of PCs as value
print "Indexing chr${chr} covariate data...\n";
my %hash=();
while (<$covs>){
	chomp;
	my @line=split(/\t/, $_);
	my $key=join("\t", @line[0 .. 1]);
	my $pcs=join("\t", @line[2 .. $#line]);
	
	$hash{$key}=$pcs;
}

# Create hash keyed by singleton positions, with input line as value
print "Writing output...\n";
my @POS;
my %poshash=();
while (<$positions>) {
	chomp;
	my @line=split(/\t/, $_);
	my $chr=$line[0];
	my $pos=$line[1];
	my $base = substr($seq, $pos-1, 1);
	
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
	
	my $key2=join("\t", $chr, $bin);
	
	# write line if site has non-N context and exists in covariate data
	if (($sequence !~ /N/) & (defined($hash{$key2}))) {
		print "getting data for position: $pos...\n";
		my $key2=join("\t", $chr, $bin);
		print "getting conservation score for: $pos...\n";
		my $score=&getCons($pos);
		print OUT "$chr\t$bin\t$pos\t$sequence\t 0 \t$hash{$key2}\t$score\n";	 
	} 
}

sub getCons{

	my $pos=shift;
	my $score;
	my $start;

	my $startpos=0;
	my $indexpos1=0;
	my $indexpos2;
	my $row=0;
	my @positions;
	while (<$cons>){
		chomp;
		# my $key;
		if($_ =~ /fixedStep/){
			my @head = split(/[=,\s]+/, $_);
			$startpos=$head[4]; # start position of next rows
			
			$row=$.; # row where header located
			$indexpos1=$pos-$startpos+1; # number of rows past header where score for site should be located, if within range
			$indexpos2=0;

		} else {
			$indexpos2++;
		}
		
		if($indexpos1==$indexpos2){
			print "Using start position: $startpos at row: $row\n";
			$score=$_;
			last;
		}
	}
	


	return $score;
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
