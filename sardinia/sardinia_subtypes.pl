#!/usr/local/bin/perl

##############################################################################
# Script outputs bcftools summary file containing genotype information
# used in process_subtypes.R to obtain subtype-specific counts
# 
# Input consists of genome-wide vcfs produced by sardinia_pipe.pl
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

my $parentdir="/net/bipolar/jedidiah/sardinia/mask";
my $bcftools="/net/bipolar/jedidiah/bcftools/bcftools";

my @files = <$parentdir/*.vcf.gz>;
foreach my $file (@files) {
	my $filename=fileparse($file);
	my $path=dirname($file);
	my $prefix = substr($filename, 0, index($filename, '.'));
	my $cmd="$bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AC[\t%SAMPLE=%GT]\n' $file > /net/bipolar/jedidiah/$prefix.summary &";
	
	# print "$cmd\n";
	
	&forkExecWait($cmd);
}

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