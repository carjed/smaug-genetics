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

my $seq="123456789012345678901234567890";

# my $mask_file = "/net/bipolar/jedidiah/reference_data/20140520.strict_mask.bed";
my $mask_file = "/net/bipolar/jedidiah/smaug-genetics/smaug-sandbox/test_mask.bed";
open my $mask, '<', $mask_file or die "can't open $mask_file: $!";

readline($mask);

while (<$mask>){
	my @line=split(/\t/, $_);
	my $c1=$line[1]-1;
	my $c2=$line[2];
	my $newseq=substr($seq, $c1, $c2-$c1, "N" x ($c2-$c1));
}

# my $c1=4;
# my $c2=8;
# my $range=$c2-$c1;
# my $rep="N" x $range;

# my $newseq=substr($seq, $c1, $c2-$c1, "N" x ($c2-$c1));

print "$seq\n";