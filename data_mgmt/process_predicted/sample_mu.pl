#!/usr/local/bin/perl

##############################################################################
# Step 0:
# Sample sites for non-mutated background
##############################################################################
use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Path qw(make_path);
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));

my $config = LoadFile("$configpath/_config.yaml");

# my $adj = $config->{adj};
my $adj = 3;
my $data = $config->{data};
my $parentdir = $config->{parentdir};
my $seed = $config->{rseed};
my $subseq = $adj*2+1;

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef getMotif);

my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );

print "Preparing de novo data...\n";
my $prepdnmcmd = "Rscript $parentdir/smaug-genetics/R/read_dnms.r TRUE $parentdir";
forkExecWait($prepdnmcmd);

srand($seed);

my $outfile = "$parentdir/output/predicted/validation_sites.txt";
open my $outFH, '>', $outfile or die "can't write to $outfile: $!\n";

foreach my $chr (1..22){

  print "Getting reference for chr$chr...\n";
  my $f_fasta;
  if($data eq "mask"){
    $f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
  } else {
    $f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
  }

  my $seq=getRef($f_fasta, $chr);
  my $seqlength=length($seq);

  foreach my $categ (@categs){
    my $predfile = "$parentdir/output/predicted/chr$chr.$categ.txt";
    open my $inFH, '<', $predfile or die "can't open $predfile: $!";

    my $nsites=`wc -l $predfile | cut -f1 -d' '`;
    chomp $nsites;

    print "Sampling chr$chr $categ sites...\n";
    while(<$inFH>){
      if(0.005>rand){
        chomp;
    		my @line=split(/\t/, $_);
        my $pos = $line[1];
    		my $localseq = substr($seq, $pos-$adj-1, $subseq);
        my $seqp = getMotif($localseq, $adj);
        print $outFH "$_\t0\t$categ\t$seqp\tall\n";
      }
    }

    my $dnmfile = "$parentdir/reference_data/DNMs/GoNL_$categ.txt";
    open my $dnmFH, '<', $dnmfile or die "can't open $dnmfile: $!";

    my $tmpfile_s = "$parentdir/reference_data/DNMs/tmp_s.txt";
    my $buildquerycmd = "grep \"\\s$chr\\s\" $dnmfile | cut -f 3 > $tmpfile_s";
    forkExecWait($buildquerycmd);

    my $tmpfile_o = "$parentdir/reference_data/DNMs/tmp_o.txt";
    my $dnmannocmd = "grep -Fwf  $tmpfile_s $predfile > $tmpfile_o";
    forkExecWait($dnmannocmd);

    print "Appending chr$chr $categ de novos...\n";
    while(<$dnmFH>){
      chomp;
      my @line=split(/\t/, $_);
      my $dnmid = $line[0];
      my $dnmchr = $line[1];
      my $dnmpos = $line[2];
      next unless $dnmchr == $chr;

      my $dnmlocalseq = substr($seq, $dnmpos-$adj-1, $subseq);
      my $dnmseqp = getMotif($dnmlocalseq, $adj);

      my $rateline = `grep -m 1 -Fw $dnmpos $tmpfile_o`;
      chomp $rateline;

      if(length($rateline)>2){
        print $outFH "$rateline\t1\t$categ\t$dnmseqp\t$dnmid\n";
      }
    }
  }
}
