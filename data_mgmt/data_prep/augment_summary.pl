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
use File::Basename;
use File::Path qw(make_path);
use Benchmark;
use FindBin;
# use lib "$FindBin::Bin";
use FaSlice;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

print "Script will run with the following parameters:\n";
for (sort keys %{$config}) {
    say "$_: $config->{$_}";
}

# my $adj = $config->{adj};
# my $adj=1;
my $mac = $config->{mac};
my $binw = $config->{binw};
my $data = $config->{data};
# my $bin_scheme = $config->{bin_scheme};
# my $bin_scheme = "band";
my $parentdir = $config->{parentdir};
my $count_motifs = $config->{count_motifs};
my $expand_summ = $config->{expand_summ};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef getMotif);

my $chr=$ARGV[0];
my $adj=$ARGV[1];

##############################################################################
#Process inputs
##############################################################################
my $subseq=2;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

my $bw=$binw/1000;

##############################################################################
# Read in files and initialize outputs
##############################################################################
# my $in_path = "/net/bipolar/jedidiah/testpipe/summaries";
# my $out_path = "$parentdir/output/${subseq}bp_${bw}k_${mac}_${data}2";
my $out_path = "$parentdir/motif_counts/$subseq-mers/$data";
make_path("$out_path");

##############################################################################
# Counts possible mutable sites per bin for 6 main categories
# and for local sequence analysis if selected
##############################################################################
if($count_motifs eq "TRUE"){
	my $start_time=new Benchmark;
	print "Counting motifs...\n";
	# print "seqlength: $length\n";
  my $fname;
  if($data eq "full"){
    $fname = "$parentdir/reference_data/human_g1k_v37/chr$chr.fasta.gz";
  } elsif($data eq "mask"){
    $fname = "$parentdir/reference_data/human_g1k_v37_mask/chr$chr.fasta.gz";
  }

  # if ( -e "$fname$chr.fasta.gz" ) { $fname = "$fname$chr.fasta.gz"; }
  my $fa;

  my $startpos;
  my $endpos;
	my @motifs;
  my $header;
  my $bin_out;

  my @schemes = qw( fixed band );

  foreach my $bin_scheme (@schemes){
    if($bin_scheme eq "fixed"){
      print "Getting fixed bins\n";
      my $fixedfile = "$parentdir/reference_data/genome.${bw}kb.sorted.bed";
      open my $fixedFH, '<', $fixedfile or die "$fixedfile: $!";
      $fa = FaSlice->new(file=>$fname, oob=>'N', size=>$binw);
      $bin_out = "$out_path/chr$chr.$subseq-mer_motifs_${data}.txt";

      $header = "CHR\tMotif\tnMotifs\n";

      readWindows($fixedFH, $bin_out, $header, $fa);

    } elsif(($bin_scheme eq "band") && ($adj==1)) {
      print "getting bands\n";
      my $bandfile = "$parentdir/reference_data/cytoBand.txt";
      open my $bandFH, '<', $bandfile or die "$bandfile: $!";
      $fa = FaSlice->new(file=>$fname, oob=>'N',size=>1_000_000);
      $bin_out = "$out_path/chr$chr.$subseq-mer_motifs_band_${data}.txt";

      $header = "CHR\tSTART\tEND\tBAND\tgieStain\tBIN\tMotif\tnMotifs\n";

      readWindows2($bandFH, $bin_out, $header, $fa);

    }
  }

	my $end_time=new Benchmark;
	my $difference = timediff($end_time, $start_time);
	print "Done. ";
	print "Runtime: ", timestr($difference), "\n";
}

sub readWindows {
  my $windowFH = shift;
  my $bin_out = shift;
  my $header = shift;
  my $fa = shift;

  my $startpos;
  my $endpos;
  my @motifs;
  my $bandno = 1;
  open(my $outFH, '>', $bin_out) or die "can't write to $bin_out: $!\n";

  print $outFH $header;

  my %full_count=();
  # my $tmpseq = $fa->get_slice($chr, 1, 1000000);
  # @motifs = ($tmpseq =~ /(?=([ACGT]{$subseq}))/g);
  @motifs = glob "{A,C,G,T}"x $subseq;
  $full_count{$_}++ for @motifs;

  while(<$windowFH>){
    chomp;
    my @line=split(/\t/, $_);
    my $chrind=$line[0];

    if($chrind eq "chr$chr"){
      $startpos = $line[1]+1;
      $endpos = $line[2]+2;

      my $binseq = $fa->get_slice($chr, $startpos, $endpos);

      @motifs = ($binseq =~ /(?=([ACGT]{$subseq}))/g);

        # get overall load
        my %tri_count=();
        $tri_count{$_}++ for @motifs;

        foreach my $motif (sort keys %tri_count) {
          my $altmotif = $motif;
          $altmotif =~ tr/ACGT/TGCA/;
          $altmotif = reverse $altmotif;

          # my $seqp = "$motif\($altmotif\)";

          my $sum;
          if(exists($tri_count{$motif}) && exists($tri_count{$altmotif})){
            $sum=$tri_count{$motif}+$tri_count{$altmotif};
          } elsif(exists($tri_count{$motif}) && !exists($tri_count{$altmotif})) {
            $sum=$tri_count{$motif};
          } elsif(!exists($tri_count{$motif}) && exists($tri_count{$altmotif})) {
            $sum=$tri_count{$altmotif};
          }

          $full_count{$motif}=$full_count{$motif}+$sum;

          # print $outFH "$first\t$bin\t$seqp\t$sum\n";
        }

      # writeCounts($_, $bandno, \@motifs, $outFH);
      $bandno = $bandno+1;
    }
  }

  foreach my $motif (sort keys %full_count){
    my $altmotif = $motif;
    $altmotif =~ tr/ACGT/TGCA/;
    $altmotif = reverse $altmotif;

    my $seqp = "$motif\($altmotif\)";
    my $sum = $full_count{$motif};
    # print $outFH "$first\t$bin\t$seqp\t$sum\n";
    print $outFH "$chr\t$seqp\t$sum\n";
  }

}

sub readWindows2 {
  my $windowFH = shift;
  my $bin_out = shift;
  my $header = shift;
  my $fa = shift;

  my $startpos;
  my $endpos;
  my @motifs;
  my $bandno = 1;
  open(my $outFH, '>', $bin_out) or die "can't write to $bin_out: $!\n";

  print $outFH $header;
  while(<$windowFH>){
    chomp;
    my @line=split(/\t/, $_);
    my $chrind=$line[0];

    if($chrind eq "chr$chr"){
      $startpos = $line[1]+1;
      $endpos = $line[2]+2;

      my $binseq = $fa->get_slice($chr, $startpos, $endpos);

      @motifs = ($binseq =~ /(?=([ACGT]{$subseq}))/g);

      writeCounts($_, $bandno, \@motifs, $outFH);
      $bandno = $bandno+1;
    }
  }
}

##############################################################################
# Read motif counts from hash table, sum counts symmetric motifs and write out
# counting strategy modified from https://www.biostars.org/p/5143/
##############################################################################
sub writeCounts {
  my $first = $_[0];
	my $bin = $_[1];
	my @motifs = @{$_[2]};
  my $outFH = $_[3];
	my %tri_count=();
	$tri_count{$_}++ for @motifs;

	foreach my $motif (sort keys %tri_count) {
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

    print $outFH "$first\t$bin\t$seqp\t$sum\n";
	}
}
