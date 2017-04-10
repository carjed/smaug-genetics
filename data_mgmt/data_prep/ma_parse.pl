#!/usr/local/bin/perl

##############################################################################
# This script reads a vcf and extracts singletons that pass specified filter(s)
# Singletons within multiallelic sites are included
#
# Under development:
#   -read both compressed and uncompressed vcfs
#   -allow multiple filters to be specified on cmd line
##############################################################################

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use Compress::Zlib;
use IO::Compress::Gzip;

my $help=0;
my $man=0;
my $invcf='';
my $outvcf='';

GetOptions ('i=s'=> \$invcf,
'o=s' => \$outvcf,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

# Read uncompressed vcf
# open my $vcf, '<', $invcf or die "can't open $invcf: $!";

# Read compressed vcf
my $vcf = gzopen($invcf, "rb") or
  die "can't open $invcf: $gzerrno";

# Initialize gzipped output
# my $OUT = open "| bgzip -c > $outvcf") or
#   die "Could not write to $outvcf: $!";

my $succ = open(my $OUT, '>', $outvcf);
$OUT = *STDOUT unless $succ;

while($vcf->gzreadline($_) > 0){
  chomp;
  if($_ =~ /^#/){
    print $OUT "$_\n";
  } else {
    my @line = split(/\s+/, $_);

    # Only process sites that pass filter
    # can add other conditions; e.g., QUAL/MQ filters
    if($line[6] eq "PASS"){

      # Get alternate allele(s)
      my $ALT = $line[4];

      # Check and process multiallelic sites
      if($ALT =~ /(,)/){

        # Coerce alternate allele field into array
        my @ALTS = split(',', $ALT);

        # Get string of allele counts from info field
        my ($AC) = $line[7] =~ /AC=([(\d+),]+(\d+));/;

        # Coerce allele counts to array
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
            print $OUT "$newline\n";
          }
        }
      } else {
        my ($AC) = $line[7] =~ /AC=(\d+);/;

        if($AC == 1){
          print $OUT "$_\n";
        }
      }
    }
  }
}

$vcf -> gzclose();
close $OUT if $succ;

__END__
=head1 NAME

ma_parse.pl: parse multiallelic singletons in vcf

=head1 SYNOPSIS

        ma_parse.pl [OPTIONS]
        Options:
		--i			input.vcf
		--o			output.vcf

=head1 OPTIONS

=over 8

=item B<--help>

Display this documentation

=item B<--i>

MANDATORY: specify input vcf

=item B<--o>

MANDATORY: specify output vcf

=back

=head1 DESCRIPTION

B<ma_parse.pl> annotates summary files and counts motifs in genome over fixed-width windows

=head1 AUTHOR

=over

Jedidiah Carlson E<10> Department of Bioinformatics E<10> University of Michigan

=cut
