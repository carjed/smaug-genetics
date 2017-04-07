#!/usr/local/bin/perl

##############################################################################
# This script annotates tab-delimited vcf summary files with sequence motif
# and mutation category information, and counts the number of mutable motifs
# in the reference genome
##############################################################################

##############################################################################
#Initialize inputs, options, and defaults and define errors
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
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname($relpath);

my $config = LoadFile("$configpath/_config.yaml");

print "Script will run with the following parameters:\n";
for (sort keys %{$config}) {
    say "$_: $config->{$_}";
}

my $adj = $config->{adj};
my $mac = $config->{mac};
my $binw = $config->{binw};
my $data = $config->{data};
my $bin_scheme = $config->{bin_scheme};
my $parentdir = $config->{parentdir};
my $count_motifs = $config->{count_motifs};
my $expand_summ = $config->{expand_summ};

use lib "$parentdir/lib";
use SmaugFunctions;

my $chr=$ARGV[0];

##############################################################################
#Process inputs
##############################################################################
my $subseq=2;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

my $bw=$binw/1000;

my $nextchr='';
if ($chr<22) {
	$nextchr=$chr+1;
} elsif ($chr==22) {
	$nextchr="X";
} else {
	$nextchr="Y";
}

##############################################################################
# Read in files and initialize outputs
# download hg37 from nih.gov if missing
##############################################################################
my $in_path = "/net/bipolar/jedidiah/testpipe/summaries";
my $out_path = "$parentdir/output/${subseq}bp_${bw}k_${mac}_${data}";
make_path("$out_path");

if($data eq "mask"){
  $f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
} else {
  $f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
}

my $seq=getRef($f_fasta, $chr);

my $seqlength=length($seq);
print "seqlength: $seqlength\n";
print "Done\n";

##############################################################################
# Counts possible mutable sites per bin for 6 main categories
# and for local sequence analysis if selected
# -also returns GC content per bin
##############################################################################
if($count_motifs eq "TRUE"){
	my $start_time=new Benchmark;
	print "Counting motifs...\n";
	# print "seqlength: $length\n";

	my $bin_flag;
	my $bin_out;
	if($bin_scheme eq "fixed"){
		$bin_out = "$out_path/chr$chr.motif_counts_fixed.txt";
	} elsif($bin_scheme eq "band") {
		$bin_out = "$out_path/chr$chr.motif_counts_band.txt";
	} else {
		$bin_out = "$out_path/chr$chr.motif_counts_all.txt";
	}

	open(BIN, '>', $bin_out) or die "can't write to $bin_out: $!\n";

	&countMotifs($bin_scheme);

	my $end_time=new Benchmark;
	my $difference = timediff($end_time, $start_time);
	print "Done. ";
	print "Runtime: ", timestr($difference), "\n";
}

##############################################################################
# Output expanded summary file based on selected options
# -passed to R script along with bins file(s)
##############################################################################
if($expand_summ eq "TRUE"){
	print "Expanding summary file...\n";
	my $start_time=new Benchmark;

	my $f_summ = "$in_path/${mac}_full/chr$chr.summary";
	open my $summ, '<', $f_summ or die "can't open $f_summ: $!";

	my $outfile = "$out_path/chr$chr.expanded.summary";
	open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

	if ($mac eq "singletons") {
		print OUT "CHR\tPOS\tREF\tALT\tDP\tAN\tSEQ\tALTSEQ\tSequence\tCategory\tCategory2\n";
	} elsif ($mac eq "doubletons") {
		print OUT "CHR\tPOS\tREF\tALT\tANNO\t";
	}

	#readline($summ); #<-throws out summary header if it exists

	while (<$summ>) {
		chomp;
		next unless $_ =~ /^[^,]*$/;

		# chomp $row;
		my @line=split(/\t/, $_);
		my $pos=$line[1];
		my $REF=$line[2];
		my $ALT=$line[3];
		my $localseq = substr($seq, $pos-$adj-1, $subseq);
		my $altlocalseq = reverse $localseq;
		$altlocalseq  =~ tr/ACGT/TGCA/;

		# keep only sites in fully parameterized motif
		if($localseq =~ /^[ACGT]+$/){

			my $ref1 = substr($localseq, $adj, 1);
			my $ref2 = substr($altlocalseq, $adj, 1);

			my $seqp;

			if($ref1 ~~ [qw( A C )]){
				$seqp = "$localseq\($altlocalseq\)";
			} else {
				$seqp = "$altlocalseq\($localseq\)";
			}

			my $CAT = "${REF}${ALT}";
			my $Category;
			if($CAT ~~ [qw( AC TG )]){
				$Category = "AT_CG";
			} elsif($CAT ~~ [qw( AG TC )]){
				$Category = "AT_GC";
			} elsif($CAT ~~ [qw( AT TA )]){
				$Category = "AT_TA";
			} elsif($CAT ~~ [qw( GA CT )]){
				$Category = "GC_AT";
			} elsif($CAT ~~ [qw( GC CG )]){
				$Category = "GC_CG";
			} elsif($CAT ~~ [qw( GT CA )]){
				$Category = "GC_TA";
			}

			my $Category2;
			if(substr($seqp, $adj, 2) eq "CG"){
				$Category2 = "cpg_$Category";
			} else {
				$Category2 = $Category;
			}
			# print OUT "$_\t$localseq\t$altlocalseq\t$gcprop\n";
			print OUT "$_\t$localseq\t$altlocalseq\t$seqp\t$Category\t$Category2\n";
		}
	}

	my $end_time=new Benchmark;
	my $difference = timediff($end_time, $start_time);
	print "Done. ";
	print "Runtime: ", timestr($difference), "\n";
}

# ##############################################################################
# # fork-exec-wait subroutine
# ##############################################################################
# sub forkExecWait {
#     my $cmd = shift;
#     #print "forkExecWait(): $cmd\n";
#     my $kidpid;
#     if ( !defined($kidpid = fork()) ) {
# 		die "Cannot fork: $!";
#     }
#     elsif ( $kidpid == 0 ) {
# 		exec($cmd);
# 		die "Cannot exec $cmd: $!";
#     }
#     else {
# 		waitpid($kidpid,0);
#     }
# }
#
# ##############################################################################
# # Get reference
# ##############################################################################
# sub getRef{
# 	my $f_fasta;
# 	if($data eq "mask"){
# 		$f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
# 	} else {
# 		$f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
# 	}
#
# 	if (-e $f_fasta) {
# 		print "Using reference genome: $f_fasta\n";
# 	} else {
# 		print "Reference genome not found in directory. Would you like to download one? (y/n): ";
# 		my $choice = <>;
# 		chomp $choice;
# 		if ($choice eq "y") {
# 			my $dlcmd="wget -P $parentdir/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz";
# 			&forkExecWait($dlcmd);
# 			my $unzipcmd="gunzip $parentdir/human_g1k_v37.fasta";
# 			&forkExecWait($unzipcmd);
# 		} else {
# 			die "Please upload an appropriate reference genome to $parentdir/reference_data/ \n";
# 		}
# 	}
#
# 	open my $fasta, '<', $f_fasta or die "can't open $f_fasta: $!";
# 	print "Getting reference sequence for chromosome $chr...\n";
#
# 	my $seq;
# 	while (<$fasta>) {
# 		chomp;
# 		if (/>$chr /../>$nextchr /) {
# 			next if />$chr / || />$nextchr /;
# 			$seq .=$_;
# 		}
# 	}
#
# 	return $seq;
# }

##############################################################################
# Subroutine counts occurrence of motifs per bin
##############################################################################
sub countMotifs{
	my $bin_flag = shift;
	my $length=length($seq);
	my $numbins=ceil($length/$binw);
	my $bin;

	print BIN "CHR\tBIN\tMOTIF\tCOUNT\n";
	my @motifs;
	if($bin_flag eq "fixed"){
		for my $i (0 .. $numbins-1) {
			@motifs=(substr($seq, $i*$binw, $binw)=~/(?=([ACGT]{$subseq}))/g);
			&writeCounts($i, \@motifs);
		}
	} elsif($bin_flag eq "band") {
		my $bandfile = "$parentdir/reference_data/cytoBand.txt";
		open my $bandFH, '<', $bandfile or die "$bandfile: $!";

		my @motifs;
		while(<$bandFH>){
			chomp;
			my @line=split(/\t/, $_);
			my $chrind=$line[0];
			if($chrind eq "chr$chr"){
				my $start=$line[1];
				my $end=$line[2];
				my $length=$end-$start;
				my $band=$line[3];
				my $stain=$line[4];

				@motifs=(substr($seq, $start, $length)=~/(?=([ACGT]{$subseq}))/g);

				### code below equivalent to writeCounts() sub
				my %tri_count=();
				# @tri_count{@a}=@b;
				$tri_count{$_}++ for @motifs;

				foreach my $motif (sort keys %tri_count) {
					# if ($count !~ /N|^G|^T/) {
					my $altmotif = $motif;
					$altmotif =~ tr/ACGT/TGCA/;
					$altmotif = reverse $altmotif;

					my $seqp = "$motif\($altmotif\)";

					my $sum;
					if(exists($tri_count{$motif}) && exists($tri_count{$altmotif})){
						$sum=$tri_count{$motif}+$tri_count{$altmotif};
					} elsif(exists($tri_count{$motif}) && !exists($tri_count{$altmotif})) {
						$sum=$tri_count{$motif};
					} elsif(!exists($tri_count{$motif}) && exists($tri_count{$altmotif})) {
						$sum=$tri_count{$altmotif};
					}

					print BIN "$_\t$seqp\t$sum\n";
					# }
				}
				###
			}
		}
	}	else {
		@motifs=($seq=~/(?=([ACGT]{$subseq}))/g);
		$bin = 0;
		&writeCounts($bin, \@motifs);
	}

	# Modified counting strategy from https://www.biostars.org/p/5143/
	# my @motifs=($seq=~/(?=([ACGT]{$subseq}))/g);
	#for(@trinucs){print "$_\n"};

	print "Done\n";
}

##############################################################################
# Read motif counts from hash table, sum counts symmetric motifs and write out
##############################################################################
sub writeCounts{
	my $bin = $_[0]+1;
	my @motifs = @{$_[1]};
	my %tri_count=();
	# @tri_count{@a}=@b;
	$tri_count{$_}++ for @motifs;

	foreach my $motif (sort keys %tri_count) {
		# if ($count !~ /N|^G|^T/) {
		my $altmotif = $motif;
		$altmotif =~ tr/ACGT/TGCA/;
		$altmotif = reverse $altmotif;

		my $seqp = "$motif\($altmotif\)";

		my $sum;
		if(exists($tri_count{$motif}) && exists($tri_count{$altmotif})){
			$sum=$tri_count{$motif}+$tri_count{$altmotif};
		} elsif(exists($tri_count{$motif}) && !exists($tri_count{$altmotif})) {
			$sum=$tri_count{$motif};
		} elsif(!exists($tri_count{$motif}) && exists($tri_count{$altmotif})) {
			$sum=$tri_count{$altmotif};
		}

		print BIN "$chr\t$bin\t$seqp\t$sum\n";
		# }
	}
}
