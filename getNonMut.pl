#!/usr/local/bin/perl

##############################################################################
# Used to obtain full data from logistic regression model
# loops through reference genome and outputs 1 line per base, as long as
# valid covariate data exists
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
use Cwd;
use Benchmark;
use Tie::File;

# Set options and inputs
my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";

my $baseopt;
my $chr;
my $categ;
my $bw = 100;
my $adj=1;

GetOptions ('b=s'=> \$baseopt,
'chr=s'=> \$chr,
'categ=s' => \$categ,
'bw=i' => \$bw,
'adj=i' => \$adj) or pod2usage(1);

my $f_covs = "$parentdir/output/logmod_data/${bw}kb_mut_cov2.txt";

make_path("$parentdir/output/logmod_data/chr${chr}/");

my $b1;
my $b2;
if($baseopt eq "AT"){
	$b1="A";
	$b2="T";
} else {
	$b1="C";
	$b2="G";
}

my $mask_flag=0;
my $binwidth=$bw*1000;

my $nextchr;
if ($chr<22) {
	$nextchr=$chr+1;
} elsif ($chr==22) {
	$nextchr="X";
}

my $subseq=2;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

# initialize output file
my $outfile = "$parentdir/output/logmod_data/chr${chr}_${categ}_sites.txt";
my $OUT;
open($OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

# initialize covariate data
open my $covs, '<', $f_covs or die "can't open $f_covs: $!";

# initialize singleton file
my $f_positions = "$parentdir/output/logmod_data/chr${chr}_${categ}_pos_examples.txt"; #main line for full processing
open my $positions, '<', $f_positions or die "can't open $f_positions: $!";

# Index motif file names
my $f_mlist = "$parentdir/output/7bp_1000k_rates.txt";
open my $mlist, '<', $f_mlist or die "can't open $f_mlist: $!";

our %fhash=();
my @fn;
while(<$mlist>){
  chomp;
  my @line=split(/\t/, $_);
  my $seq=$line[1];
  my $cat=$line[2];

  if($cat eq $categ){
    #my $key=join("\t", @line[0 .. 1]);
    #my $pcs=join("\t", @line[2 .. $#line]);

    my $filename="$parentdir/output/logmod_data/chr${chr}/chr${chr}_${categ}_$seq.txt";
    push(@fn, $filename);
    $fhash{$seq}=$filename;
    # print "$hash{$_}\n";
  }
}

my %handles = get_write_handles(@fn);

# initialize phastCons data
# my $f_cons = "$parentdir/reference_data/chr$chr.phastCons46way.primates.wigFix";
# open my $cons, '<', $f_cons or die "can't open $f_cons: $!";

# Get reference sequence
my $seq=&getRef();
my $altseq=$seq;
$altseq =~ tr/ACGT/TGCA/;

my $seqlength=length($seq);
# print "seqlength of chr$chr: $max\n"; #<-used to validate that getRef() returns correct seq length

my $printheader=0;
if($printheader==1){
	print $OUT "CHR \t POS \t BIN \t Sequence \t mut \n"; #<-add header to output, if needed
}

# Create hash keyed by Chr/Bin pairs, with row of PCs as value
print "Indexing chr${chr} covariate data...\n";
our %hash=();
while (<$covs>){
	chomp;
	my @line=split(/\t/, $_);
	my $key=join("\t", @line[0 .. 1]);
	my $pcs=join("\t", @line[2 .. $#line]);

	$hash{$key}=$pcs;
}

# my $key=join("\t", 20, 100);
# print "$hash{$key}\n";

# Create hash keyed by singleton positions, with input line as value
print "Indexing chr${chr}: ${categ} singleton positions...\n";
my @POS;
my %poshash=();
while (<$positions>) {
	chomp;
	my @line=split(/\t/, $_);
	my $key=$line[2];
	push (@POS, $key);

	$poshash{$key}=$_;
}

# Write data files
print "Writing chr${chr}: ${categ} data file...\n";
for my $strpos (0 .. $seqlength){
	my $base = substr($seq, $strpos, 1);
	my $pos = $strpos+1;

	my $bin = ceil($pos/$binwidth);
	my $key2=join("\t", $chr, $bin);

	if(defined $hash{$key2}){
		if(($base =~ /$b1|$b2/) & (!exists $poshash{$pos})){
			# push (@POS, $pos); # add position to exclusion list
			my $localseq = substr($seq, $pos-$adj-1, $subseq);
			my $altlocalseq = reverse substr($altseq, $pos-$adj-1, $subseq);

			# Coerce local sequence info to format used in R
			my $sequence;
			if(substr($localseq,$adj,1) lt substr($altlocalseq,$adj,1)){
				$sequence = $localseq . '(' . $altlocalseq . ')';
			} else {
				$sequence = $altlocalseq . '(' . $localseq . ')';
			}

			# write line if site has non-N context
			if ($sequence !~ /[MNSW]/) {
				my $covs=&updateCovs($chr, $bin, $pos);

				my $file=$fhash{$sequence};
				# my $mref={$handles{$file}};

				print $OUT "$chr\t$bin\t$pos\t$sequence\t 0 \t$covs\n";
				print {$handles{$file}} "$chr\t$bin\t$pos\t$sequence\t 0 \t$covs\n";
			}
		}elsif(exists $poshash{$pos}){
			my $covs=&updateCovs($chr, $bin, $pos);
			my @line=split(/\t/, $poshash{$pos});
			my $sequence=$line[3];

			my $file=$fhash{$sequence};
			# my $mref={$handles{$file}};

			print $OUT "$poshash{$pos}\t$covs\n";
			print {$handles{$file}} "$poshash{$pos}\t$covs\n";
		}
	}
}

print "Done\n";

sub getRef{
	my $f_fasta;
	if($mask_flag==1){
		$f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
	} else {
		$f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
	}

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

	##############################################################################
	# Retrieve reference sequence for selected chromosome
	# -also returns symmetric sequence to be used in local sequence analysis
	##############################################################################

	print "Getting reference sequence for chromosome $chr...\n";

	my $seq;
	if($mask_flag==1){
		while (<$fasta>) {
			chomp;
			if (/^>$chr$/../^>$nextchr$/) {
				next if /^>$chr$/ || /^>$nextchr$/;
				$seq .=$_;
			}
		}
	} else {
		while (<$fasta>) {
			chomp;
			if (/>$chr /../>$nextchr /) {
				next if />$chr / || />$nextchr /;
				$seq .=$_;
			}
		}
	}

	return $seq;
}


sub updateCovs{

	my $CHR=shift;
	my $BIN=shift;
	my $pos=shift;

	# Get keys for current and neighboring bins
	my $c_linekey=join("\t", $CHR, $BIN);
	my $p_linekey=join("\t", $CHR, $BIN-1);
	my $n_linekey=join("\t", $CHR, $BIN+1);

	# my $pos=$line[2];
	my $posmin=$BIN*$binwidth-$binwidth;

	# Get relative position
	my $relpos = ($pos-$posmin)/$binwidth;

	# Calculate proportion
	my $prop_c_bin = -abs($relpos-0.5)+1;

	# Get covariates of containing bin
	my $o_line = $hash{$c_linekey};

	# print "$c_linekey\n";
	# print "$o_line\n";

	my @c_feats = split(/\t/, $o_line);

	# Calculate covariates in current bin proportional to position
	foreach my $x (@c_feats) { $x = $x * $prop_c_bin; }

	my @sum;
	# Repeat for adjacent window
	if(($relpos-0.5<=0) && exists($hash{$p_linekey})){
		my $prop_p_bin = -$relpos+0.5;
		my @p_feats = split(/\t/, $hash{$p_linekey});
		foreach my $x (@p_feats) { $x = $x * $prop_p_bin; }
		@sum = pairwise { $a + $b } @c_feats, @p_feats;

	} elsif(($relpos-0.5>0) && exists($hash{$n_linekey})){
		my $prop_n_bin = $relpos-0.5;
		my @n_feats = split(/\t/, $hash{$n_linekey});
		foreach my $x (@n_feats) { $x = $x * $prop_n_bin; }
		@sum = pairwise { $a + $b } @c_feats, @n_feats;
	} else{
		@sum = @c_feats;
	}

	my $covs=join("\t", @sum[0 .. $#sum]);
	# print "$covs\n";
	return $covs;
}

# Subroutine reads array of filenames and returns file handles
sub get_write_handles {
  my @file_names = @_;
  my %file_handles;
  foreach (@file_names) {
    open my $fh, '>', $_ or next;
    $file_handles{$_} = $fh;
  }
  return %file_handles;
}
