#!/usr/local/bin/perl

##############################################################################
# Script takes list of chr/positions, annotates with rate and k-mer according
# to specified file
#
# Currently used to annotate rocdat_comb.txt as rocdat_comb_3bp.txt, which is
# passed to the roc.r script
##############################################################################

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Benchmark;
use Math::Round;

my $parentdir="/net/bipolar/jedidiah/mutation";

my $help=0;
my $man=0;
# my $chr;
my $adj=1;
# my $f_rates;
my $f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
my $f_positions;
my $f_rates;
my $seqflag;

GetOptions (
# 'chr=i'=> \$chr,
'adj=i' => \$adj,
# 'f_rates=s' => \$f_rates,
'ref=s' => \$f_fasta,
'in=s' => \$f_positions,
'seq' => \$seqflag,
'rates=s' => \$f_rates,
'help|?'=> \$help,
'ref=s' => \$f_fasta,
man => \$man) or pod2usage(1);

my $subseq = $adj*2+1;

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;

my %cathash = (
    AT_CG => 0,
    AT_GC => 1,
    AT_TA => 2,
		GC_AT => 3,
    GC_CG => 4,
    GC_TA => 5,
);

# Initialize and hash rate table
# my $f_rates = "$parentdir/ERV_${subseq}bp_rates.txt";
open my $rates, '<', $f_rates or die "can't open $f_rates: $!";
readline($rates); #<-throws out header
# print "Hashing rate table...\n";
my %hash=();
while (<$rates>){
	chomp;
	my @line=split(/\t/, $_);
	my $key=$line[0];
	# my $vals=join("\t", nearest(0.0001, @line[1 .. $#line]));
	# my $vals=sum(@line[1 .. $#line]);
	$hash{$key}=[@line[1 .. $#line]];
}

# my $f_positions="/net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_full_mask.txt";
# my $f_positions="$parentdir/output/predicted/full/rocdat_comb_sort.txt";
# my $f_positions="$parentdir/output/predicted/full/rocdat_comb_7bp.txt";
# my $f_positions="$parentdir/output/rocdat.txt";
# my $f_positions="/net/bipolar/jedidiah/mutation/reference_data/uk10k_dnms_s.txt";
open my $positions, '<', $f_positions or die "can't open $f_positions: $!";

# Initialize output file
# my $outfile ="$parentdir/output/predicted/full/rocdat_comb2_${subseq}bp.txt";
my $outfile ="$parentdir/output/rocdat.${subseq}bp.txt";
# my $outfile ="$parentdir/reference_data/uk10kdnms_3mer.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";
# print OUT "CHR\tPOS\tAT_CG\tAT_GC\tAT_TA\tGC_AT\tGC_CG\tGC_TA\n";

# Get reference sequence
# my $seq=&getRef();
# my $altseq=$seq;
# $altseq =~ tr/ACGT/TGCA/;

# print "Indexing sites in chr${chr}...\n";
my @sites;
my $prevchr=0;
my $seq;
my $altseq;
readline($positions);
while(<$positions>){
	chomp;
	my $linestr=$_;
	my @line=split(/\t/, $_);
	my $sitechr=$line[0];
	my $pos=$line[1];
	my $categ=$line[5];
	my $catind=$cathash{$categ};

	if($sitechr!=$prevchr){
		$seq=&getRef($sitechr);
		$altseq=$seq;
		$altseq =~ tr/ACGT/TGCA/;
	}

	my $base=substr($seq, $pos, 1);

	my $localseq = substr($seq, $pos-$adj-1, $subseq);

	# print "$localseq\n";

	if($localseq!~/[MNSW]/){
		my $altlocalseq = reverse substr($altseq, $pos-$adj-1, $subseq);

		my $sequence;
		if(substr($localseq,$adj,1) lt substr($altlocalseq,$adj,1)){
			$sequence = $localseq . '(' . $altlocalseq . ')';
		} else {
			$sequence = $altlocalseq . '(' . $localseq . ')';
		}

		# print OUT "$chr\t$i\t$hash{$sequence}\n";
		if($seq){
			print OUT "$linestr\t$sequence\t$hash{$sequence}[$catind]\t\n";
		}	else {
			print OUT "$linestr\t$hash{$sequence}[$catind]\t\n";
		}
	}

	$prevchr=$sitechr;
}


sub getRef{
	my $chr=shift;
	my $nextchr;
	if ($chr<22) {
		$nextchr=$chr+1;
	} elsif ($chr==22) {
		$nextchr="X";
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

__END__
=head1 NAME

anno_rate.pl - Whole-genome annotation utility for 5bp motif relative mutation rates

=head1 SYNOPSIS

        anno_rate.pl [OPTIONS]
        Options:
		--help			program documentation
		--f-rates		/path/to/5bp_rates.txt
		--ref			/path/to/human_g1k_v37.fasta


=head1 OPTIONS

=over 8

=item B<--help>

Display this documentation

=item B<--f-rates>

Path to rate table. Default is same directory as this script.

=item B<--ref>

Path to reference genome.  Default is same directory as this script. Tested with ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz

=back

=head1 DESCRIPTION

B<anno_rate.pl> Whole-genome annotation utility for 5bp motif relative mutation rates.

=head1 AUTHOR

=over

Jedidiah Carlson E<10> Department of Bioinformatics E<10> University of Michigan

=cut
