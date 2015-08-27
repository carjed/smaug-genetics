#!/usr/local/bin/perl

##############################################################################
# Annotate summary file with CpG status
# used to exclude these sites from trio analysis
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

my $parentdir="/net/bipolar/jedidiah";

##############################################################################
#Process inputs
##############################################################################

my $cov;
GetOptions ('cov=i'=> \$cov);

my $subseq=2;

##############################################################################
# Read in files and initialize outputs
# download hg37 from nih.gov if missing
# -Will eventually update summary file location to match pipe.pl
##############################################################################

my $f_fasta = "$parentdir/human_g1k_v37.fasta";

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

my $f_summ = "/net/bipolar/jedidiah/mask2/${cov}x_full_mask.summary";

open my $summ, '<', $f_summ or die "can't open $f_summ: $!";

my $outfile = "/net/bipolar/jedidiah/mask2/expanded/${cov}x_expanded.summary";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";


print OUT "CHR\tPOS\tREF\tALT\tDP\tAC\tS28009\tS28038\tS28047\tPAIR\tGC\n";


##############################################################################
# Output expanded summary file based on selected options
# -passed to R script along with bins file(s)
##############################################################################

print "Creating data file...\n";

my @POS;
my @NEWSUMM;

while (<$summ>) {
	push (@POS, (split(/\t/, $_))[1]);
	push (@NEWSUMM, $_);
}

my $chr=1;
my $nextchr=2;
my $seq;

while (<$fasta>) {
	chomp;
	if (/>$chr/../>$nextchr/) {
		next if />$chr/ || />$nextchr/;
		$seq .=$_;
	}
}

foreach my $row (@NEWSUMM) {
	chomp $row;
	my @line=split(/\t/, $row);
	my $pos=$line[1];
	
	my $prevchr=$chr;
	$chr=$line[0];
	
	if ($chr<22) {
		$nextchr=$chr+1;
	} elsif ($chr==22) {
		$nextchr="X";
	}
	
	if ($chr!=$prevchr){
		my $seq='';
		while (<$fasta>) {
			chomp;
			if (/>$chr/../>$nextchr/) {
				next if />$chr/ || />$nextchr/;
				$seq .=$_;
			}
		}
	}
		
	
		
	my $localseq = substr($seq, $pos-1, $subseq);
	# my $altseq=$localseq;
	# $altseq =~ tr/ACGT/TGCA/;
	
	my $nextseq = substr($seq, $pos, $subseq);
	# my $altnextseq=$nextseq;
	# $altnextseq =~ tr/ACGT/TGCA/;
	
	my $cpg;
	if($localseq eq "CG" or $localseq eq "GC" or $nextseq eq "CG" or $nextseq eq "GC"){
		$cpg=1;
	}
	else{
		$cpg=0;
	}
	
	# $seq='';
	my $gcprop = &getGC($pos);

	# print OUT "$row\t$localseq\t$gcprop\n";
	print OUT "$row\t$cpg\t$gcprop\n";
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
# GC Content Subroutine
# Returns proportion of G or C nucleotides in reference genome for given
# flanking region of site
##############################################################################
sub getGC {
	my $site=shift;
	# my $binseq=shift;
	my $binwidth=100;
	my $region_s=$site-($binwidth/2);
	#my $region_e=$site-$binwidth/2;
	my $region=substr($seq,$region_s,$binwidth);
	
	my $abase=($region =~ tr/A//);
	my $cbase=($region =~ tr/C//);
	my $gbase=($region =~ tr/G//);
	my $tbase=($region =~ tr/T//);

	my $gcsum=$cbase+$gbase;
	my $total=$abase+$cbase+$gbase+$tbase;
	my $gc_content=0.5;
	
	if ($total != 0) {
		$gc_content = $gcsum/$total;
	}
	return $gc_content;
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
		--cpg			CpG site analysis?
		--hot			recombination hotspots?
		--anno			annotation(s)

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

=item B<--cpg>

toggles extra analysis specific to CpG sites

=item B<--hot>

toggles extra analysis for distance to nearest recombination hotspot 

=item B<--anno>

comma-separated list consisting of any of the following annotations:

Intergenic
Intron
Nonsynonymous
Exon
Synonymous
Utr3
Utr5
Upstream
Downstream
Stop_Gain
Stop_Loss
Start_Loss
Essential_Splice_Site
Normal_Splice_Site

=back

=head1 DESCRIPTION

B<my-prog.pl> is doing something.

=head1 AUTHOR

=over

Jedidiah Carlson E<10> Department of Biostatistics E<10> University of Michigan

=cut