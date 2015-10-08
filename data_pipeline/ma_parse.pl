#!/usr/local/bin/perl

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
use Compress::Zlib;

my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";

my $help=0;
my $man=0;
my $invcf='';
my $outvcf='';

GetOptions ('i=s'=> \$invcf,
'o=s' => \$outvcf,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

# open my $vcf, '<', $invcf or die "can't open $invcf: $!";
my $vcf = gzopen($invcf, "rb") or die "can't open $invcf: $gzerrno";

open(OUT, '>', $outvcf) or die "can't write to $outvcf: $!\n";

while($vcf->gzreadline($_) > 0){
  chomp;
  if($_ =~ /^#/){
    print OUT "$_\n";
  } else {
    my @line = split(/\s+/, $_);
    my $ALT = $line[4];

    # Check and process multiallelic sites
    if($ALT =~ /(,)/){

      # Get coerce alternate allele field into array
      my @ALTS = split(',', $ALT);

      # Get string of allele counts from info field
      my ($AC) = $line[7] =~ /AC=([(\d+),]+(\d+));/;

      # coerce allele counts to array
      my @ACS = split(',', $AC);

      # Loop through alts and extract singletons
      for my $i (0 .. $#ACS){
        my $iAC = $ACS[$i];

        if($iAC == 1){
          my $oldline = $_;
          my $maAC = "AC=$AC;";
          my $newAC = "AC=1;";
          my $newALT = $ALTS[$i];

          (my $newline = $oldline) =~ s/$maAC/$newAC/;
          $newline =~ s/$ALT/$newALT/g;
          print OUT "$newline\n";
        }
      }
    } else {
      my ($AC) = $line[7] =~ /AC=(\d+);/;

      if($AC == 1){
        print OUT "$_\n";
      }
    }
  }
}

$vcf -> gzclose();

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

=item B<--data>

specify whether to use summaries from all singletons (full) or those that pass the strict filters (strict)

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

B<ref5.pl> annotates summary files and counts motifs in genome over fixed-width windows

=head1 AUTHOR

=over

Jedidiah Carlson E<10> Department of Bioinformatics E<10> University of Michigan

=cut
