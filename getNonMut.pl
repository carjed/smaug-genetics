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
my $parentdir="/net/bipolar/jedidiah/mutation";

my $baseopt;
my $chr;
my $categ;

GetOptions ('b=s'=> \$baseopt,
'chr=s'=> \$chr,
'categ=s' => \$categ) or pod2usage(1);

my $b1;
my $b2;
if($baseopt eq "AT"){
	$b1="A";
	$b2="T";
} else {
	$b1="C";
	$b2="G";
}

my $mask_flag=1;
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

my $outfile = "/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/chr${chr}_${categ}_negative_examples.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

my $f_covs = "/net/bipolar/jedidiah/mutation/smaug-genetics/sandbox/mut_cov2.txt"; #main line for full processing
open my $covs, '<', $f_covs or die "can't open $f_covs: $!";

# Create hash keyed by Chr/Bin pairs, with row of PCs as value
my %hash=();
while (<$covs>){
	my @line=split(/\t/, $_);
	my $key=join("\t", @line[0 .. 1]);
	my $pcs=join("\t", @line[2 .. $#line]);
	
	# print "$key \t $pcs \n";
	
	$hash{$key}=$pcs;
	
	# print "$hash{$key}\n";
}

# my $testkey="20\t153";
# print "$hash{$testkey}\n";

# my $seq="ACGTAAGCCGTAAC";
# 2,3,7,8,9,10,14

my $seq=&getRef();
my $altseq=$seq;
$altseq =~ tr/ACGT/TGCA/;

my $min=5;
my $max=length($seq)-5;

print "seqlength of chr$chr: $max\n";

# print OUT "CHR \t POS \t BIN \t Sequence \t mut \n";

# Read in file with positions to exclude from analyses
my @POS;

my $f_positions = "/net/bipolar/jedidiah/mutation/smaug-genetics/chr${chr}_${categ}_exclusion_list.txt"; #main line for full processing
open my $positions, '<', $f_positions or die "can't open $f_positions: $!";

print "Reading chr${chr}: ${categ} exclusion list...\n";

while (<$positions>) {
	# next unless $_ =~ /^[^,]*$/;
	push (@POS, $_);
}

print "Done\n";

print "Writing chr${chr}: ${categ} data file...\n";
for my $i ($min .. $max){
# for my $i (1 .. 100) {

	# my $pos = $min + int(rand($max - $min));
	# my $base = substr($seq, $pos-1, 1);
	
	my $pos = $i;
	my $base = substr($seq, $pos-1, 1);

	if(($base =~ /$b1|$b2/) & ($pos !~ @POS)){
		# push (@POS, $pos); # add position to exclusion list
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
		
		if (($sequence !~ /N/) & (defined($hash{$key2}))) {
			my $key2=join("\t", $chr, $bin);
			# my $pcs2=$hash{$key2};
			
			# print "Key: $key2\n";
			# print "$pcs2\n";
		
			print OUT "$chr\t$bin\t$pos\t$sequence\t 0 \t $hash{$key2}";	
		} 
	}
}
print "Done\n";

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