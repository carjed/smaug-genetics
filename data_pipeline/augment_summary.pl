#!/usr/local/bin/perl

##############################################################################
# SMAUG: Singleton Mutation Analysis Utility with Graphics
##############################################################################
# SMAUG uses extremely rare variants to visualize changes in mutation rates
# across the genome. The ref5.pl script takes bcftools summary files and
#
# Jedidiah Carlson
# Department of Bioinformatics
# The University of Michigan
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

my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";

my $help=0;
my $man=0;
my $chr;
my $mac;
my $binwidth=100000;
my $adj=0;
my $data="full";
my $mask_flag='';
my $count_motifs='';
my $by_bin='';
my $expand_summ='';

GetOptions ('chr=i'=> \$chr,
'mac=i'=> \$mac,
'b=i' => \$binwidth,
'adj=i' => \$adj,
'data=s' => \$data,
'mf' => \$mask_flag,
'motifs' => \$count_motifs,
'bybin' => \$by_bin,
'summ' => \$expand_summ,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;

if (!$chr | !$mac) {
	pod2usage("$0: Missing mandatory argument.");
}

##############################################################################
#Process inputs
##############################################################################
my $macl;
if ($mac==1) {
	$macl = "singletons";
}
if ($mac==2) {
	$macl = "doubletons";
}
if ($mac==11){
	$macl="common";
}

my $subseq=2;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

my $bw=$binwidth/1000;

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
my $out_path = "$parentdir/output/${subseq}bp_${bw}k_${macl}_${data}";
make_path("$out_path");

my $seq=&getRef();

my $seqlength=length($seq);
print "seqlength: $seqlength\n";
print "Done\n";

##############################################################################
# Counts possible mutable sites per bin for 6 main categories
# and for local sequence analysis if selected
# -also returns GC content per bin
##############################################################################
if($count_motifs){
	my $start_time=new Benchmark;
	print "Counting motifs...\n";
	# print "seqlength: $length\n";

	my $bin_flag;
	my $bin_out;
	if($by_bin){
		$bin_flag = "bybin";
		$bin_out = "$out_path/chr$chr.motif_counts_binned.txt";
	} else {
		$bin_flag = "bychr";
		$bin_out = "$out_path/chr$chr.motif_counts.txt";
	}

	open(BIN, '>', $bin_out) or die "can't write to $bin_out: $!\n";
	&countMotifs($bin_flag);
	# &binCounts();

	my $end_time=new Benchmark;
	my $difference = timediff($end_time, $start_time);
	print "Done. ";
	print "Runtime: ", timestr($difference), "\n";
}

##############################################################################
# Output expanded summary file based on selected options
# -passed to R script along with bins file(s)
##############################################################################
if($expand_summ){
	print "Expanding summary file...\n";
	my $start_time=new Benchmark;

	my $f_summ = "$in_path/${macl}_${data}/chr$chr.summary";
	open my $summ, '<', $f_summ or die "can't open $f_summ: $!";

	my $outfile = "$out_path/chr$chr.expanded.summary";
	open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

	if ($mac==1) {
		print OUT "CHR\tPOS\tREF\tALT\tDP\tAN\tSEQ\tALTSEQ\tSequence\tCategory\tCategory2\n";
	} elsif ($mac==2) {
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

##############################################################################
# fork-exec-wait subroutine
##############################################################################
sub forkExecWait {
    my $cmd = shift;
    #print "forkExecWait(): $cmd\n";
    my $kidpid;
    if ( !defined($kidpid = fork()) ) {
		die "Cannot fork: $!";
    }
    elsif ( $kidpid == 0 ) {
		exec($cmd);
		die "Cannot exec $cmd: $!";
    }
    else {
		waitpid($kidpid,0);
    }
}

##############################################################################
# Get reference
##############################################################################
sub getRef{
	my $f_fasta;
	if($mask_flag){
		$f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
	} else {
		$f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
	}

	if (-e $f_fasta) {
		print "Using reference genome: $f_fasta\n";
	} else {
		print "Reference genome not found in directory. Would you like to download one? (y/n): ";
		my $choice = <>;
		chomp $choice;
		if ($choice eq "y") {
			my $dlcmd="wget -P $parentdir/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz";
			&forkExecWait($dlcmd);
			my $unzipcmd="gunzip $parentdir/human_g1k_v37.fasta";
			&forkExecWait($unzipcmd);
		} else {
			die "Please upload an appropriate reference genome to $parentdir/reference_data/ \n";
		}
	}

	open my $fasta, '<', $f_fasta or die "can't open $f_fasta: $!";
	print "Getting reference sequence for chromosome $chr...\n";

	my $seq;
	while (<$fasta>) {
		chomp;
		if (/>$chr /../>$nextchr /) {
			next if />$chr / || />$nextchr /;
			$seq .=$_;
		}
	}

	return $seq;
}

##############################################################################
# Subroutine counts occurrence of motifs per bin
##############################################################################
sub countMotifs{
	my $bin_flag = shift;
	my $length=length($seq);
	my $numbins=ceil($length/$binwidth);
	my $bin;

	print BIN "CHR\tBIN\tMOTIF\tCOUNT\n";
	my @motifs;
	if($bin_flag eq "bybin"){
		for my $i (0 .. $numbins-1) {
			@motifs=(substr($seq, $i*$binwidth, $binwidth)=~/(?=([ACGT]{$subseq}))/g);
			&writeCounts($i, \@motifs);
		}
	} else {
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
	my $bin = $_[0];
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

__END__
=head1 NAME

ref5.pl - SMAUG: Singleton Mutation Analysis Utility with Graphics

=head1 SYNOPSIS

        ref5.pl [OPTIONS]
        Options:
		--help			program documentation
		--chr			chromosome
		--mac			minor allele count
		--b			binwidth
		--adj			number of adjacent nucleotides
		--data			data subset to use

=head1 OPTIONS

=over 8

=item B<--help>

Display this documentation

=item B<--chr>

MANDATORY: specify chromosome for analysis

=item B<--mac>

MANDATORY: specify minor allele count of sites in existing summary file

=item B<--b>

specify bin width for histograms (default is 100,000)

=item B<--adj>

specify number of adjacent nucleotides in either direction from the variant to include in analysis
default includes only the adjacent 3' nucleotide for CpG distinction

=item B<--data>

specify whether to use summaries from all singletons (full) or those that pass the strict filters (strict)

=back

=head1 DESCRIPTION

B<ref5.pl> annotates summary files and counts motifs in genome over fixed-width windows

=head1 AUTHOR

=over

Jedidiah Carlson E<10> Department of Bioinformatics E<10> University of Michigan

=cut
