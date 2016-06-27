#!/usr/local/bin/perl

##############################################################################
# Scans summarized glf files (every 10bp; colnames: chr, pos, ref, dp) and
# outputs mean depth
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
use Math::Round;
use Cwd;
use Benchmark;
use Tie::File;

# Set options and inputs
my $help=0;
my $man=0;
my $index;

GetOptions ('ind=i'=> \$index,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;

my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";
my $glfdir="$parentdir/output/glf_depth/chr$chr/";

my $filelist="/net/bipolar/jedidiah/mutation/output/glf_depth/glf_filelist.txt";
open my $files, '<', $filelist or die "$filelist: $!";

my $indfile;
while( <$files> ) {
    if($. == $index) {
        $indfile = $_;
        last;
    }
}

# open my $ind, '<', $indfile or die "$indfile: $!";

my @filepath=split m%/%, $indfile;

my $fname="$filepath[7]/$filepath[8]/$filepath[9].dp";

print "$_\n" foreach(@filepath);

my $glfcmd="samtools-hybrid glfview $indfile | cut -f1-4 | awk '\$2%10==0 && \$3 ~ /[ACGT]/' > /net/bipolar/jedidiah/mutation/output/glf_depth/$fname";
print "Following command will be run:\n";
print "$glfcmd\n";
# &forkExecWait($glfcmd);

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
