#!/usr/local/bin/perl

##############################################################################
# SMAUG: Singleton Mutation Analysis Utility with Graphics
#
# Jedidiah Carlson
# Department of Biostatistics
# The University of Michigan
# Last Revision: 06/10/2014
#
##############################################################################
# SUMMARY:
# SMAUG uses extremely rare variants to visualize changes in mutation rates
# across the genome. The following genomic features are currently implemented
# or are in development:
# -local sequence (+/- 1 adjacent nucleotide)
# -CpG status
# -Distance to nearest recombination hotspot
# -Fuctional annotation
# -GC content
# -Average read depth
# 
# CURRENT RESTRICTIONS:
# -Only works with singleton summary files (--mac 1); expanded doubleton
# 	summary files are not yet created
# -only takes tri-nucleotide sequences (--adj 1)
#
# TO-DO:
# -update directories here and in R script to be consistent with pipe.pl
# -Integrate with pipe.pl
# -stitch together heatmaps for multiple chromosomes
# -subset annotations
# -directly integrate call to bcftools for summary files?
# -create directory for per-chromosome sequence files
#
# LONG-TERM GOALS:
# -Cross-reference with annotation data
# -Integrate somatic mutation info
# -Better integration of different frequency class info (e.g. shared plots)
# -Hidden Markov Model
#
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

my $wdir=getcwd;
my $parentdir=dirname($wdir);

my $help=0;
my $man=0;
my $chr;
my $mac;
my $binwidth=100000;
my $adj=0;
my $cpg='';
my $hot='';
my $anno='';

GetOptions ('chr=i'=> \$chr,
'mac=i'=> \$mac,
'b=i' => \$binwidth,
'adj=i' => \$adj,
'cpg' => \$cpg,
'hot' => \$hot,
'anno=s' => \$anno,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;

if (!$chr | !$mac) {
	pod2usage("$0: Missing mandatory argument.");
}

print "Local subsequence and CpG command entered simultaneously--overriding CpG analysis\n" if ($cpg && $adj!=0);

##############################################################################
#Process mandatory inputs
##############################################################################

make_path("$parentdir/images/chr$chr");
my $imgdir="$parentdir/images/chr$chr";

my @annos=qw(Intergenic Intron Nonsynonymous Exon Synonymous Utr3 Utr5 Upstream Downstream Stop_Gain Stop_Loss Start_Loss Essential_Splice_Site Normal_Splice_Site);

if ($anno && $anno ~~ @annos) {
	make_path("$parentdir/images/chr$chr/anno/$anno");
	$imgdir="$parentdir/images/chr$chr/anno/$anno";
} elsif ($anno && !($anno~~@annos)) {
	print "Invalid annotation entered--all sites will be used.  See help file for further information.\n";
}

print "$imgdir\n";

my $macl;
if ($mac==1) {
	$macl = "singletons";
}
if ($mac==2) {
	$macl = "doubletons";
}

my $cpg_flag;
if ($cpg && $adj==0) {
	$cpg_flag="on";
} else {
	$cpg_flag="off";
}

my $hot_flag;
if ($hot) {
	$hot_flag="on";
} else {
	$hot_flag="off";
}

my $subseq=2;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

my $nextchr;
if ($chr<22) {
	$nextchr=$chr+1;
} elsif ($chr==22) {
	$nextchr="X";
}

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

#my $f_summ = "/net/bipolar/jedidiah/bcftools/summaries/$macl/all/chr$chr.$macl.summary.txt";
my $f_summ = "$parentdir/testpipe/summaries/chr$chr.summary";
open my $summ, '<', $f_summ or die "can't open $f_summ: $!";

my $outfile = "expanded.summary";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

my $bin_out = 'bin_out.txt';
open(BIN, '>', $bin_out) or die "can't write to $bin_out: $!\n";

my $bin_out2 = 'bin_out2.txt';
open(BIN2, '>', $bin_out2) or die "can't write to $bin_out2: $!\n";

##############################################################################
# Retrieve reference sequence for selected chromosome
# -also returns symmetric sequence for local sequence analysis
##############################################################################

print "Getting reference sequence for chromosome $chr...\n";

my $seq;
while (<$fasta>) {
	chomp;
	if (/>$chr/../>$nextchr/) {
		next if />$chr/ || />$nextchr/;
		$seq .=$_;
	}
}

my $altseq=$seq;
$altseq =~ tr/ACGT/TGCA/;

# Validate overall GC content of chromosome
my $abase=($seq =~ tr/A//);
my $cbase=($seq =~ tr/C//);
my $gbase=($seq =~ tr/G//);
my $tbase=($seq =~ tr/T//);

my $gcsum=$cbase+$gbase;
my $total=$abase+$cbase+$gbase+$tbase;
my $gc = $gcsum/$total;

print "GC content of chr$chr: $gc\n";

print "Done\n";

##############################################################################
#call CpGI script
##############################################################################

my @cpgi_index;
my $temp_fasta;

if ($cpg && $adj==0) {
	$temp_fasta = 'temp.fasta';
	open(TEMP, '>', $temp_fasta) or die "can't write to $temp_fasta: $!\n";
	print TEMP ">chr$chr\n";
	print TEMP "$seq\n";
	
	print "Running CpG Island analysis...\n";
	my $cpgicmd="perl CpGcluster.pl temp.fasta 50 1E-5 > CpGCluster.log";
	&forkExecWait($cpgicmd);
	print "Done\n";
	
	my $f_cpgi = "$wdir/temp.cpg";
	open my $cpgi, '<', $f_cpgi or die "can't open $f_cpgi: $!";

	while (<$cpgi>) {
		push(@cpgi_index, $_);
	}
}

##############################################################################
# Counts possible mutable sites per bin for 6 main categories
# and for local sequence analysis if selected
# -also returns GC content per bin, to be implemented later
##############################################################################

print "Getting bin counts...\n";

if ($adj==1) {

	my $length=length($seq);
	my $numbins=ceil($length/$binwidth);
	my $bin;
	my @A;
	my @C;
	my @G;
	my @T;
	my @a= glob "{A,C,G,T}"x $subseq;
	my @b = (0) x (scalar @a);

	my @trinucs=($seq=~/(?=(.{$subseq}))/g);
	my %tri_count=();
	@tri_count{@a}=@b;
	$tri_count{$_}++ for @trinucs;

	print BIN2 "SEQ\tCOUNT\n";
	foreach my $count (sort keys %tri_count) {
		if ($count !~ /N/) {
			print BIN2 "$count\t$tri_count{$count}\n";
		} 
	}

	print BIN "AT\tCG\tprop_GC\tBIN\t";
	foreach my $tri (@a) {
		my $alttri=$tri;
		$alttri =~ tr/ACGT/TGCA/;
		my $min = minstr($tri, $alttri);
		my $max = maxstr($tri, $alttri);
		if ($tri !~ /^G|^T/) {
			print BIN "$min($max)\t";
		}
	}
	print BIN "\n";

	for my $i (0 .. $numbins-1) {
		$A[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/A//);
		$C[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/C//);
		$G[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/G//);
		$T[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/T//);

		my @trinucs=(substr($seq, $i*$binwidth, $binwidth)=~/(?=(.{$subseq}))/g);
		#for(@trinucs){print "$_\n"};
		my %tri_count=();
		@tri_count{@a}=@b;
		$tri_count{$_}++ for @trinucs;
		
		my $sum_at=$A[$i]+$T[$i];
		my $sum_cg=$C[$i]+$G[$i];
		my $GC=0;
		if (($sum_at+$sum_cg)!=0) {
			$GC=($sum_cg/($sum_at+$sum_cg));
		}
		
		$bin=$i+1;

		print BIN "$sum_at\t$sum_cg\t$GC\t$bin\t";
		#print BIN "$_:$tri_count{$_}\t" for sort keys(%tri_count);
		foreach my $count (sort keys %tri_count) {
			if ($count !~ /N|^G|^T/) {
				my $altcount= $count;
				$altcount =~ tr/ACGT/TGCA/;
				my $sum=$tri_count{$count}+$tri_count{$altcount};
				print BIN "$sum\t";
			} 
		}
		print BIN "\n";
	}
} else {
	my $length=length($seq);
	my $numbins=ceil($length/$binwidth);
	my $bin;
	my @A;
	my @C;
	my @G;
	my @T;

	print BIN "AT\tCG\tprop_GC\tBIN\n";

	for my $i (0 .. $numbins-1) {
		$A[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/A//);
		$C[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/C//);
		$G[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/G//);
		$T[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/T//);

        my $sum_at=$A[$i]+$T[$i];
        my $sum_cg=$C[$i]+$G[$i];
        my $GC=0;
        
		if (($sum_at+$sum_cg)!=0) {
            $GC=($sum_cg/($sum_at+$sum_cg));
        }
        
        $bin=$i+1;

        print BIN "$sum_at\t$sum_cg\t$GC\t$bin\n";
	}
}

print "Done\n";

##############################################################################
# Get hotspot coordinates if selected
##############################################################################

my @loci;
if ($hot) {
	my $f_hotspots = "$parentdir/genetic_map/hotspots.txt";
	#my $f_hotspots="test_hotspots.txt";
	open my $hotspots, '<', $f_hotspots or die "can't open $f_hotspots: $!";

	print "Getting hotspot coordinates...\n";
	
	readline($hotspots); #<-throws out header
	while (<$hotspots>) {
		if ($_ =~ /^chr$chr/) {
			push (@loci, $_);
		}
	}
	print "Done\n";
}

##############################################################################
# Output expanded summary file based on selected options
# -passed to R script along with bins file(s)
##############################################################################

print "Creating data file...\n";

my @POS;
my @NEWSUMM;
#readline($summ); #<-throws out summary header if it exists

if ($anno && $anno~~@annos) {
	print "subsetting $anno sites\n";
	while (<$summ>) {
		push (@POS, (split(/\t/, $_))[1]);
		if ($_ =~ /$anno/) {
			push (@NEWSUMM, $_);
		}
	}
} else {
	while (<$summ>) {
		push (@POS, (split(/\t/, $_))[1]);
		push (@NEWSUMM, $_);
	}
}

if ($cpg && $adj==0) {
	
	print OUT "CHR\tPOS\tREF\tALT\tDP\tAN\tANNO\tPAIR\tCPGI\n";

	foreach my $row (@NEWSUMM) {
		chomp $row;
		my @line=split(/\t/, $row);
		my $pos=$line[1];
		my $localseq = substr($seq, $pos-$adj-1, $subseq);
		my $hit=0;	
		foreach my $cpgi_int (@cpgi_index) {
			my @pair = split(/,/, $cpgi_int);
			my $start=$pair[0];
			my $end=$pair[1];
			
			if (($pos >= $start) && ($pos <= $end)) {
				$hit=1;
				last;
			}
		}

		print OUT "$row\t$localseq\t$hit\n";
	}
} elsif ($hot) {

	print OUT "CHR\tPOS\tREF\tALT\tDP\tAN\tANNO\tSEQ\tALTSEQ\tDIST\n";

	foreach my $row (@NEWSUMM) {
		chomp $row;
		my @line=split(/\t/, $row);
		my $pos=$line[1];
		my $distance = &dist2Hotspot($pos);
		my $localseq = substr($seq, $pos-$adj-1, $subseq);
		my $altlocalseq = substr($altseq, $pos-$adj-1, $subseq);
		
		print OUT "$row\t$localseq\t$altlocalseq\t$distance\n";
	}
} else {

	print OUT "CHR\tPOS\tREF\tALT\tDP\tAN\tANNO\tSEQ\tALTSEQ\n";

	foreach my $row (@NEWSUMM) {
		chomp $row;
		my @line=split(/\t/, $row);
		my $pos=$line[1];
		my $localseq = substr($seq, $pos-$adj-1, $subseq);
		my $altlocalseq = substr($altseq, $pos-$adj-1, $subseq);
		
		print OUT "$row\t$localseq\t$altlocalseq\n";
	}
}

print "Done\n";

##############################################################################
#Run selected R script
##############################################################################

my $args="$chr $macl $binwidth $cpg_flag $outfile $adj $hot_flag $imgdir";
my $cmd="Rscript prop.R $args";

print "Running R script...\n";
&forkExecWait($cmd);
print "Done. See images folder for output.\n";

##############################################################################
#Clean up temp files
##############################################################################

if ($cpg && $adj==0) {
	unlink $temp_fasta;
}

unlink $outfile;
unlink $bin_out;
unlink $bin_out2;

my $plots_out="Rplots.pdf";
unlink $plots_out;

my $Rlog="R.log";
unlink $Rlog;

my $CpGlog="CpGCluster.log";
unlink $CpGlog;

my $CpGlog2="temp.cpg-log.txt";
unlink $CpGlog2;

my $CpGtemp="temp.cpg";
unlink $CpGtemp;

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
# Hotspot Subroutine
# Input: variant position
# Output: Distance to hotspot
##############################################################################

sub dist2Hotspot {
	my $site = shift;
	my $prevabsstart=10e10;
	my $prevabsend=10e10;
	my $prevdist=10e10;
	
	my $dist;
	my $distout;
	
	#returns distance to first locus if site occurs before 
	#returns distance to last locus if site occurs after
	my @first = split(/\t/, $loci[0]);
	my @last = split(/\t/, $loci[$#loci]);
	
	if ($site < $first[2]) {
		$dist = $first[2]-$site;
		return $dist;
		last;
	} elsif ($site > $last[3]) {
		$dist = $site-$last[3];
		return $dist;
		last;	
	}
	
	foreach my $locus (@loci) {
		my @elements=split(/\t/, $locus);
		my $absstart=abs($site-$elements[2]);
		my $absend=abs($site-$elements[3]);

		# returns 0 and stops if site is in hotspot
		# checks distance to current hotspot; if previous hotspot is closer $prevdist is returned
		# if site is near last hotspot in list (with $dist<$prevdist), loop ends and $dist is returned

		if ($site >=$elements[2] && $site <= $elements[3]) {
			$dist=0;
			return $dist;
			last;
		} elsif ($site < $last[3]) {
			$dist = min $absstart, $absend, $prevabsstart, $prevabsend;

			if ($prevdist <= $dist) {
				$dist=$prevdist;
				return $prevdist;
				last;
			} 
		}
		$prevabsstart=$absstart;
		$prevabsend=$absend;
		$prevdist=$dist;		
	}
	return $dist;
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
		--anno			annotation

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

subset summary file by one of the following annotations and perform analyses:

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