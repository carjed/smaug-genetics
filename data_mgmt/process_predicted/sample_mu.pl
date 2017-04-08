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

my $adj = $config->{adj};
my $mac = $config->{mac};
my $binw = $config->{binw};
my $data = $config->{data};
my $bin_scheme = $config->{bin_scheme};
my $parentdir = $config->{parentdir};
my $count_motifs = $config->{count_motifs};
my $expand_summ = $config->{expand_summ};

my $subseq=$adj*2+1;

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef getMotif);

# my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );
my @categs = qw(AT_CG);

print "Preparing de novo data...\n";
my $prepdnmcmd = "Rscript $parentdir/smaug-genetics/R/read_dnms.r TRUE $parentdir";
forkExecWait($prepdnmcmd);

# print if 0.025 > rand while <$input>;
srand(36087318);

foreach my $categ (@categs){
  my $outfile = "$parentdir/output/predicted/${categ}.sub_new2.txt";
  open my $outFH, '>>', $outfile or die "can't write to $outfile: $!\n";

  # foreach my $chr (1..22){
  foreach my $chr (21..22){

    print "Getting reference for chr$chr...\n";
    my $f_fasta;
    if($data eq "mask"){
      $f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
    } else {
      $f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
    }

    my $seq=getRef($f_fasta, $chr);
    my $predfile = "$parentdir/output/predicted/chr$chr.$categ.txt";
  	open my $inFH, '<', $predfile or die "can't open $predfile: $!";

    while(<$inFH>){
      if(0.05>rand){
        chomp;
    		my @line=split(/\t/, $_);
        my $pos = $line[1];
    		my $localseq = substr($seq, $pos-$adj-1, $subseq);
        my $seqp = getMotif($localseq, $adj);
        print $outFH "$_\t0\t$categ\t$seqp\n";
      }
    }

    my $dnmfile = "$parentdir/reference_data/DNMs/GoNL_$categ.txt";
    open my $dnmFH, '<', $predfile or die "can't open $predfile: $!";

    while(<$dnmFH>){
      chomp;
      my @line=split(/\t/, $_);
      my $dnmchr = $line[1];
      my $dnmpos = $line[2];
      next unless $dnmchr =~ $chr;

      my $localseq = substr($seq, $dnmpos-$adj-1, $subseq);
      my $seqp = getMotif($localseq, $adj);
      print "DNM line: $_\n";
      my $rateline = `grep -w $dnmpos $predfile`;
      chomp $rateline;
      print "grep result: $rateline\n";
      my @ratelinearr = split(/\t/, $rateline);
      if($dnmpos~~$rateline){
        print $outFH "$rateline\t1\t$categ\t$seqp\n";
      }
    }

    # my $tmpfile = "$parentdir/reference_data/DNMs/tmp.txt";
    #
    # my $buildquerycmd = "grep \"\\s$i\\s\" $dnmfile | cut -f 3 > $tmpfile";
    # forkExecWait($buildquerycmd);
    #
    # my $dnmannocmd = "grep -Fwf  $tmpfile $parentdir/output/predicted/chr$i.${categ}.txt | awk -v categ=\"$categ\" '{print \$0\"\\t\"1\"\\t\"categ}' >> $parentdir/reference_data/DNMs/GoNL_${categ}.anno.txt";
    # forkExecWait($dnmannocmd);

  }

}



# foreach my $categ (@categs) {
#   # Get random selection of sites with predicted rates
#   # Must modify to include category
#   print "Subsetting $categ sites...\n";
#   my $samplecmd = "awk -v categ=\"$categ\" 'BEGIN {srand()} !/^\$/ { if (rand() <= .005) print \$0\"\\t\"0\"\\t\"categ}' $parentdir/output/predicted/chr*.${categ}.txt > $parentdir/output/predicted/${categ}.sub_new.txt";
#   forkExecWait($samplecmd);
#
#   # remove per-category DNM data if it already exists
#   my $cleanupcmd = "rm -f $parentdir/reference_data/DNMs/GoNL_${categ}.anno.txt";
#   forkExecWait($cleanupcmd);
#
#   # Testing: get de novo data in same script
#   for my $i (1 .. 22) {
#     my $tmpfile = "$parentdir/reference_data/DNMs/tmp.txt";
#
#     my $buildquerycmd = "grep \"\\s$i\\s\" $parentdir/reference_data/DNMs/GoNL_${categ}.txt | cut -f 3 > $tmpfile";
#     forkExecWait($buildquerycmd);
#
#     my $dnmannocmd = "grep -Fwf  $tmpfile $parentdir/output/predicted/chr$i.${categ}.txt | awk -v categ=\"$categ\" '{print \$0\"\\t\"1\"\\t\"categ}' >> $parentdir/reference_data/DNMs/GoNL_${categ}.anno.txt";
#     forkExecWait($dnmannocmd);
#   }
# }
#
# print "Combining and sorting data...\n";
# # Combine null sites with DNMs
# my $sortcmd = "cat $parentdir/output/predicted/*.sub_new.txt $parentdir/reference_data/DNMs/*.anno.txt | sort -V -k1,2 > $parentdir/output/rocdat.sort_new.txt";
# forkExecWait($sortcmd);
