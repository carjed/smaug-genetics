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
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils 'pairwise';
use Cwd;
use Benchmark;
use Tie::File;
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

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

my $baseopt;
my $chr;
my $categ;
my $bw = 100;

GetOptions ('b=s'=> \$baseopt,
'chr=s'=> \$chr,
'categ=s' => \$categ,
'bw=i' => \$bw) or pod2usage(1);

make_path("$parentdir/output/logmod_data/chr${chr}/");

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
my $binwidth=$bw*1000;

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
my $outfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_sites.txt";

# Initialize gzipped output
open(my $OUT, "| gzip -c > $outfile") or
  die "Could not write to $outfile: $!";
# open($OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

# initialize singleton file
my $f_positions = "$parentdir/output/logmod_data/chr${chr}_${categ}_pos_examples.txt"; #main line for full processing
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
	print $OUT "CHR \t POS \t BIN \t Sequence \t mut \n"; #<-add header to output, if needed
}

# Create hash keyed by singleton positions, with input line as value
print "Indexing chr${chr}: ${categ} singleton positions...\n";
my @POS;
my %poshash=();
while (<$positions>) {
	chomp;
	my @line=split(/\t/, $_);
	my $key=$line[2];
	push (@POS, $key);

	$poshash{$key}=$_;
}

# Write data files
print "Writing chr${chr}: ${categ} data file...\n";
for my $strpos (0 .. $seqlength){
	my $base = substr($seq, $strpos, 1);
	my $pos = $strpos+1;

	my $bin = ceil($pos/$binwidth);
	my $key2=join("\t", $chr, $bin);

	# if(defined $hash{$key2}){
		if(($base =~ /$b1|$b2/) & (!exists $poshash{$pos})){
			# push (@POS, $pos); # add position to exclusion list
			my $localseq = substr($seq, $pos-$adj-1, $subseq);
			my $altlocalseq = reverse $localseq;
			$altlocalseq  =~ tr/ACGT/TGCA/;

			# Coerce local sequence info to format used in R
			my $sequence;
			if(substr($localseq,$adj,1) lt substr($altlocalseq,$adj,1)){
				$sequence = $localseq . '(' . $altlocalseq . ')';
			} else {
				$sequence = $altlocalseq . '(' . $localseq . ')';
			}

			# write line if site has non-N context
			if ($sequence =~ /\A[acgt\(\)]+\z/i) {
				print $OUT "$chr\t$bin\t$pos\t$sequence\t0\n";
			}
		} elsif(exists $poshash{$pos}){
			print $OUT "$poshash{$pos}\n";
		}
	# }
}

print "Done\n";


##############################################################################
# Subroutine reads array of filenames and returns file handles
##############################################################################
sub get_write_handles {
  my @file_names = @_;
  my %file_handles;
  foreach (@file_names) {
    open my $fh, '>', $_ or next;
    $file_handles{$_} = $fh;
  }
  return %file_handles;
}
