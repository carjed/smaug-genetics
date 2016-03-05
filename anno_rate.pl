#!/usr/local/bin/perl

##############################################################################
# Script used to annotate every non-N base in a specified chromosome with
# relative mutation rates
#
# Output is a text file with the following 8 columns:
# CHR	POS	AT_CG	AT_GC	AT_TA	GC_AT	GC_CG	GC_TA
#
# For each base position in the file, 3 rates will be output; the others will
# be NA values
#
# Currently must specify a single chromosome when running; can be modified
# to loop through all chromosomes
#
# Use perl anno_rate.pl --help for additional details
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

my $parentdir=getcwd;

my $help=0;
my $man=0;
my $chr;
my $f_rates = "$parentdir/5bp_rates.txt";
my $f_fasta = "$parentdir/human_g1k_v37.fasta";

GetOptions ('chr=i'=> \$chr,
'f_rates=s' => \$f_rates,
'ref=s' => \$f_fasta,
'help|?'=> \$help,
'ref=s' => \$f_fasta,
man => \$man) or pod2usage(1);

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;

if (!$chr) {
	pod2usage("$0: Missing mandatory argument.");
}

my $nextchr;
if ($chr<22) {
	$nextchr=$chr+1;
} elsif ($chr==22) {
	$nextchr="X";
}

# Initialize and hash rate table
open my $rates, '<', $f_rates or die "can't open $f_rates: $!";
readline($rates); #<-throws out header
print "Hashing rate table...\n";
my %hash=();
while (<$rates>){
	chomp;
	my @line=split(/\t/, $_);
	my $key=$line[0];
	# my $vals=join("\t", nearest(0.0001, @line[1 .. $#line]));
	my $vals=sum(@line[1 .. $#line]);
	$hash{$key}=$vals;
}

my $f_positions="/net/bipolar/jedidiah/mutation/output/predicted/chr${chr}_full_mask.txt";
open my $positions, '<', $f_positions or die "can't open $f_positions: $!";

my @sites;
while(<$positions>){
	chomp;
	my @line=split(/\t/, $_);
	my $pos=$line[0];
	push(@sites, $pos);
}

# Initialize output file
my $outfile ="$parentdir/chr$chr.rates.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";
# print OUT "CHR\tPOS\tAT_CG\tAT_GC\tAT_TA\tGC_AT\tGC_CG\tGC_TA\n";

# Get reference sequence
my $seq=&getRef();
my $altseq=$seq;
$altseq =~ tr/ACGT/TGCA/;

# Define motif length
my $adj=1;
my $subseq=1;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

# Get rates for each base
my $start_time=new Benchmark;
print "Writing data file...\n";
# for my $i (2 .. length($seq)-1){
foreach my $i (@sites){
	my $base=substr($seq, $i, 1);

	my $localseq = substr($seq, $i-$adj-1, $subseq);

	if($localseq!~/[MNSW]/){
		my $altlocalseq = reverse substr($altseq, $i-$adj-1, $subseq);

		my $sequence;
		if(substr($localseq,$adj,1) lt substr($altlocalseq,$adj,1)){
			$sequence = $localseq . '(' . $altlocalseq . ')';
		} else {
			$sequence = $altlocalseq . '(' . $localseq . ')';
		}

		# print OUT "$chr\t$i\t$hash{$sequence}\n";
		print OUT "$i\t$hash{$sequence}\n";
	}
}
my $end_time=new Benchmark;
my $difference = timediff($end_time, $start_time);
print "Done. ";
print "Finished in: ", timestr($difference), "\n";

sub getRef{
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
		--chr			chromosome
		--f-rates		/path/to/5bp_rates.txt
		--ref			/path/to/human_g1k_v37.fasta


=head1 OPTIONS

=over 8

=item B<--help>

Display this documentation

=item B<--chr>

MANDATORY: specify chromosome for analysis

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
