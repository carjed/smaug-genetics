#!/usr/local/bin/perl

##############################################################################
# Used to obtain full data from logistic regression model
# loops through reference genome and outputs 1 line per base, as long as
# valid covariate data exists
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

my $adj = 3;
my $data = $config->{data};
my $parentdir = $config->{parentdir};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

my $chr;
my $categ;

GetOptions ('chr=s'=> \$chr,
'categ=s' => \$categ);

make_path("$parentdir/output/logmod_data/chr${chr}/");

my $b1;
my $b2;
if(substr($categ, 0, 2) eq "AT"){
	$b1="A";
	$b2="T";
} elsif(substr($categ, 0, 2) eq "GC") {
	$b1="C";
	$b2="G";
}

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
my $outfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_full.txt.gz";

# Initialize gzipped output
open(my $OUT, "| gzip -c > $outfile") or
  die "Could not write to $outfile: $!";
# open($OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

# initialize singleton file
my $f_positions = "$parentdir/output/logmod_data/chr${chr}_${categ}_sites.txt"; #main line for full processing
open my $positions, '<', $f_positions or die "can't open $f_positions: $!";

# Index motif file names
my $f_mlist = "$parentdir/output/7bp_1000k_rates.txt";
open my $mlist, '<', $f_mlist or die "can't open $f_mlist: $!";

print "Getting reference for chr$chr...\n";
my $f_fasta;
if($data eq "mask"){
  $f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
} else {
  $f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
}

my $seq=getRef($f_fasta, $chr);
my $seqlength=length($seq);

my $printheader=0;
if($printheader==1){
	print $OUT "CHR\tPOS\tSequence\tmut\n"; #<-add header to output, if needed
}

# Create hash keyed by singleton positions, with input line as value
print "Indexing chr${chr}: ${categ} singleton positions...\n";
my @POS;
my %poshash=();
while (<$positions>) {
	chomp;
	my @line=split(/\t/, $_);
	my $key=$line[1];
	push (@POS, $key);

	$poshash{$key}=$_;
}

# hash depth file
my $dpdir="$parentdir/output/glf_depth/meandp";
my %dphash=();
my $chrfile="$dpdir/chr$chr.dp";
my $chrFH;
open($chrFH, '<', $chrfile) or
	die "Unable to open file $chrfile : $!";

while(my $dp=<$chrFH>){
	chomp($dp);
	my @dpline=split(/\t/, $dp);
	my $dp_pos=$dpline[0];
	my $depth=$dpline[1];

	$dphash{$dp_pos}=$depth;
}

# Write data files
print "Writing chr${chr}: ${categ} data file...\n";
for my $strpos ($adj .. $seqlength){
	my $base = substr($seq, $strpos, 1);
	my $pos = $strpos+1;
	my $poslim=rounddown($pos,10);

	my $outline;
	if(($base =~ /$b1|$b2/) & (!exists $poshash{$pos})){
		# push (@POS, $pos); # add position to exclusion list
		my $localseq = substr($seq, $pos-$adj-1, $subseq);
		my $altlocalseq = reverse $localseq;
		$altlocalseq  =~ tr/ACGT/TGCA/;

		# Coerce local sequence info to format used in R
		if ($localseq =~ /\A[acgt\(\)]+\z/i) {
			my $sequence;
			if(substr($localseq,$adj,1) lt substr($altlocalseq,$adj,1)){
				$sequence = $localseq . '(' . $altlocalseq . ')';
			} else {
				$sequence = $altlocalseq . '(' . $localseq . ')';
			}

		# write line if site has non-N context

			$outline = "$chr\t$pos\t$sequence\t0\t";
		}
	} elsif(exists $poshash{$pos}){
		$outline = "$poshash{$pos}\t";
	}

	# query depth hash and write if value exists
	if(exists($dphash{$poslim}) && defined($outline)){
		my $dpout=$dphash{$poslim};
		print $OUT "$outline\t$dpout\n";
	}


}

print "Done\n";

##############################################################################
# Rounding subroutines
##############################################################################
sub roundup {
  my $num = shift;
  my $roundto = shift || 1;

  return int(ceil($num/$roundto))*$roundto;
}

sub rounddown {
  my $num = shift;
  my $roundto = shift || 1;

  return int(floor($num/$roundto))*$roundto;
}
