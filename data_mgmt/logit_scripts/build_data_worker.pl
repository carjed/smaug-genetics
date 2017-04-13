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

# my $catind = $ARGV[0]-1;
#
# my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );
#
# my $categ = $categs[$catind];
#
# my $b1;
# my $b2;
# if(substr($categ, 0, 2) eq "AT"){
# 	$b1="A";
# 	$b2="T";
# } elsif(substr($categ, 0, 2) eq "GC") {
# 	$b1="C";
# 	$b2="G";
# }

# Initialize gzipped output
# my $outfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_full.txt.gz";
# open(my $outFH, "| gzip -c > $outfile") or
#   die "Could not write to $outfile: $!";

# foreach my $chr (reverse(1 .. 22)){
foreach my $chr (1 .. 22){

	# Create hash keyed by singleton positions, with input line as value
	# print "Indexing chr${chr}: ${categ} singleton positions...\n";

	# my $posfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_sites.txt";
	my $posfile = "$parentdir/output/logmod_data/chr${chr}_sites.txt";
	print "Indexing chr${chr} singleton file: $posfile...\n";

	open my $posFH, '<', $posfile or die "can't open $posfile: $!";
	my %poshash=();

	while (my $posline=<$posFH>) {
		chomp($posline);
		my @posfields=split(/\t/, $posline);
		my $key = $posfields[1];

		$poshash{$key} = $posline;
	}

	print "Getting chr${chr}: reference sequence...\n";
	my $fname;
	if($data eq "full"){
		$fname = "$parentdir/reference_data/human_g1k_v37/chr$chr.fasta.gz";
	} elsif($data eq "mask"){
		$fname = "$parentdir/reference_data/human_g1k_v37_mask/chr$chr.fasta.gz";
	}

	my $fa = FaSlice->new(file=>$fname, oob=>'N', size=>1_000_000);

	my $fixedfile = "$parentdir/reference_data/genome.5000kb.sorted.bed";
	open my $fixedFH, '<', $fixedfile or die "$fixedfile: $!";
	# $fa = FaSlice->new(file=>$fname, oob=>'N', size=>$binw);

	my $chunkpath = "$parentdir/output/logmod_data/chr$chr";
	make_path($chunkpath);

	# my $splitpath = "$parentdir/output/logmod_data/motifs3/$categ";
	my $splitpath = "$parentdir/output/logmod_data/motifs3";
	make_path($splitpath);

	my $i = 1;
	while(my $binline=<$fixedFH>){
		chomp($binline);
		my @binfields = split(/\t/, $binline);
		my $chrind = $binfields[0];
		my $startpos = $binfields[1]+1;
		my $endpos = $binfields[2];
		my $length = $endpos-$startpos+1;

		# my $outfile = "$chunkpath/$chrind.$startpos-$endpos.${categ}.txt";
		my $outfile = "$chunkpath/$chrind.$startpos-$endpos.txt";

		if($chrind eq "chr$chr"){

			# hash depth file
			my $dpdir = "$parentdir/output/glf_depth/meandp";
			# /net/bipolar/jedidiah/mutation/output/glf_depth/meandp/chr20.5000001.10000000.txt
			my $dpfile = "$dpdir/chr$chr.$startpos.$endpos.txt";
			print "Indexing chunk $i depth file: $dpfile...\n";
			open my $dpFH, '<', $dpfile or die "Unable to open file $dpfile : $!";
			my %dphash=();

			while(my $dpline=<$dpFH>){
				chomp($dpline);
				my @dpfields = split(/\t/, $dpline);
				my $dp_pos = $dpfields[0];
				my $depth = $dpfields[1];

				$dphash{$dp_pos}=$depth;
			}

			print "Writing chunk $i output file: $outfile...\n";
			open my $outFH, '>', $outfile or die "Could not write to $outfile: $!";

			my $loopstart = $startpos;
			my $loopend = $endpos;

			if($startpos <= $adj){
				$loopstart = $startpos+$adj;
			} elsif($length != 5000000) {
				$loopend = $endpos-$adj;
			}

			for my $pos ($loopstart .. $loopend){
				# my $base = $fa->get_base($chr, $pos);
				my $motif = $fa->get_slice($chr, $pos-$adj, $pos+$adj);
				my $base = substr($motif, 3, 1);

				my $outline;
				my $poslim = rounddown($pos,10);

				if(($motif =~ /\A [ACGT()]+\z/ix) && (!exists $poshash{$pos})){
					my $fullmotif = getMotif($motif, $adj);
					$outline = "$chr\t$pos\t$fullmotif\t0\t0\t0\t0\t0\t0\t";
				} elsif(exists $poshash{$pos}){
					$outline = "$poshash{$pos}\t";
				}

				# query depth hash and write if value exists
				if(exists($dphash{$poslim}) && defined($outline)){
					my $dpout = $dphash{$poslim};
					my @cols = split(/\t/, $outline);
					print $outFH "$outline\t$dpout\n";
				}
			}

			close $outFH or warn $! ? "Error closing: $!" : "Exit status $? ";

			print "Splitting chunk $i output file: $outfile...\n";
			# my $fullfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_full.txt.gz";
			# my $subcmd = "sort -k3 $outfile | awk '{print >> \"$splitpath/${categ}_\" substr(\$3, 1, 7) \".txt\"}'";
			my $subcmd = "sort -k3,3 $outfile | awk '{print >> \"$splitpath/\" substr(\$3, 1, 7) \".txt\"}'";
			# my $subcmd = "cat $outfile | awk '{print >> \"$splitpath/${categ}_\" substr(\$3, 1, 7) \".txt\"}'";
			forkExecWait($subcmd);
			$i++;
		}
	}

	print "Done\n";
}

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
