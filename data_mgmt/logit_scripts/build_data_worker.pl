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
use FaSlice;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

my $adj = 3;
my $subseq = $adj*2+1;
my $data = $config->{data};
my $parentdir = $config->{parentdir};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getMotif);

my $chr;
my $categ;

GetOptions ('chr=s'=> \$chr,
'categ=s' => \$categ);

my $b1;
my $b2;
if(substr($categ, 0, 2) eq "AT"){
	$b1="A";
	$b2="T";
} elsif(substr($categ, 0, 2) eq "GC") {
	$b1="C";
	$b2="G";
}

# Initialize gzipped output
# my $outfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_full.txt.gz";
# open(my $outFH, "| gzip -c > $outfile") or
#   die "Could not write to $outfile: $!";

# initialize singleton file
my $posfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_sites.txt";
open my $posFH, '<', $posfile or die "can't open $posfile: $!";

# Create hash keyed by singleton positions, with input line as value
print "Indexing chr${chr}: ${categ} singleton positions...\n";
my @POS;
my %poshash=();
while (<$posFH>) {
	chomp;
	my @line=split(/\t/, $_);
	my $key=$line[1];
	push (@POS, $key);

	$poshash{$key}=$_;
}

# hash depth file
my $dpdir="$parentdir/output/glf_depth/meandp";
my %dphash=();
my $dpfile="$dpdir/chr$chr.dp";
open my $dpFH, '<', $dpfile or die "Unable to open file $dpfile : $!";

while(my $dp=<$dpFH>){
	chomp($dp);
	my @dpline=split(/\t/, $dp);
	my $dp_pos=$dpline[0];
	my $depth=$dpline[1];

	$dphash{$dp_pos}=$depth;
}

print "Getting reference for chr$chr...\n";
my $fname;
if($data eq "full"){
	$fname = "$parentdir/reference_data/human_g1k_v37/chr$chr.fasta.gz";
} elsif($data eq "mask"){
	$fname = "$parentdir/reference_data/human_g1k_v37_mask/chr$chr.fasta.gz";
}

my $fa = FaSlice->new(file=>$fname, oob=>'N', size=>1_000_000);

# Write data files
print "Writing chr${chr}: ${categ} data file...\n";

my $fixedfile = "$parentdir/reference_data/genome.1000kb.sorted.bed";
open my $fixedFH, '<', $fixedfile or die "$fixedfile: $!";
# $fa = FaSlice->new(file=>$fname, oob=>'N', size=>$binw);

my $outpath = "$parentdir/output/logmod_data/motifs3/$categ";
make_path($outpath);

while(<$fixedFH>){
	chomp;
	my @line=split(/\t/, $_);
	my $chrind=$line[0];

	if($chrind eq "chr$chr"){
		my $startpos = $line[1]+1;
		my $endpos = $line[2];
		for my $pos ($startpos .. $endpos){
			my $base = $fa->get_base($chr, $pos);
			my $poslim=rounddown($pos,10);
			my $outline;

			if(($base =~ /$b1|$b2/) & (!exists $poshash{$pos})){
				my $motif = $fa->get_slice($chr, $pos-$adj, $pos+$adj);
        $motif = getMotif($motif, $adj);

				# write line if site has non-N context
				if($motif =~ /\A [ACGT()]+\z/ix){
					$outline = "$chr\t$pos\t$motif\t0\t";
				}
			} elsif(exists $poshash{$pos}){
				$outline = "$poshash{$pos}\t";
			}

			# query depth hash and write if value exists
			if(exists($dphash{$poslim}) && defined($outline)){
				my $dpout=$dphash{$poslim};

				my @cols = split(/\t/, $outline);
				my $motif = $cols[2];
				$motif = substr($motif, 0, 7);

				my $outfile = "$outpath/${categ}_$motif.txt.gz";
				open(my $outFH, "| gzip -c > $outfile") or
				  die "Could not write to $outfile: $!";

				print $outFH "$outline\t$dpout\n";
			}
		}
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
