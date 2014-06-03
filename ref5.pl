#!/usr/local/bin/perl
use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use Cwd;

##############################################################################
# SMAUG: Singleton Mutation Analysis Utility with Graphics
#
# Compares summary file from bcftools with reference genome and 
# outputs adjacent sites, then passes output to R script to plot distribution
#
# To-Do:
# -update directories here and in R scripts to be consistent with pipe.pl
# -make folder in images directory for each chromosome
# -Parallel handling of fasta start/endpoints
# -Integrate with pipe.pl
# -create top level R script with shared cmds and branch only unique cmds
# -update R scripts to account for local sequence
# -update image directories used in R scripts
# -directly integrate call to bcftools for summary files?
# -create directory for per-chromosome sequence files
#
# Long-term goals:
# -Cross-reference with annotation data and plot
# -Integrate somatic mutation info
# -Better integration of different frequency class info
# -Hidden Markov Model
#
# Associated Files (need to update for new CPGI script):
# -prop.R
# -count.R
# -bin_out.txt
# -cpg_out.txt
# -cpgi130.pl
# -human_g1k_v37.fasta
# -temp.fasta
# -
##############################################################################

##############################################################################
#Initialize inputs/options/defaults and define errors
##############################################################################
my $wdir=getcwd;
my $parentdir=dirname($wdir);

my $help=0;
my $man=0;
my $chr;
my $mac;
my $script;
my $binwidth=100000;
my $adj=0;
my $cpg='';

my $subseq=2;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

GetOptions ('chr=i'=> \$chr,
'mac=i'=> \$mac,
'script=i' => \$script,
'b=i' => \$binwidth,
'f=i' => \$adj,
'cpg' => \$cpg,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;

if (!$chr | !$mac | !$script) {
	pod2usage("$0: Missing mandatory argument.");
}

##############################################################################
#Process mandatory inputs
##############################################################################
my $macl;
if ($mac==1) {
	$macl = "singletons";
}
if ($mac==2) {
	$macl = "doubletons";
}

#my $summ = "/net/bipolar/jedidiah/bcftools/summaries/$macl/all/chr$chr.$macl.summary.txt";
my $summ = "/net/bipolar/jedidiah/testpipe/summaries/chr$chr.summary";

my $cpg_flag;
if ($cpg) {
	$cpg_flag="on";
} else {
	$cpg_flag="off";
}

my $cmd;
my $args="$chr $macl $binwidth $cpg_flag $summ";
if ($script==1) {
	$cmd="Rscript count.R $args";
}
if ($script==2) {
	$cmd="Rscript prop.R $args";
}

my $sindex=$chr-1;
my $eindex=$chr;

##############################################################################
# Read in reference fasta file and summary file and initialize outputs
# -Will eventually update summary file location to match pipe.pl
##############################################################################
my $file = "$parentdir/human_g1k_v37.fasta";

if (-e $file) {
	print "Using reference genome: $file\n";
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

open my $fasta, '<', $file or die "can't open $file: $!";
open my $index, '<', $summ or die "can't open $summ: $!";

my $cpg_out = 'cpg_out.txt';
open(OUT, '>', $cpg_out) or die "can't write to $cpg_out: $!\n";

my $bin_out = 'bin_out.txt';
open(BIN, '>', $bin_out) or die "can't write to $bin_out: $!\n";

##############################################################################
# Find all headings from fasta file to define endpoints
# Re-read in reference fasta file and subset for selected chromosome
##############################################################################
print "Getting sequence for chromosome $chr...\n";
my @endpts;
while (<$fasta>) {
	if ($_=~ />/) {
		push (@endpts, $.);
	}
}

my $start;
my $end;

for my $i (0 .. $#endpts) {
        if ($i == $sindex) {
                $start=$endpts[$i];
        }
        if ($i == $eindex) {
                $end=$endpts[$i];
        }
} 

my $file2 = "$parentdir/human_g1k_v37.fasta";
open my $fasta2, '<', $file2 or die "can't open $file: $!";

my $seq;
while (<$fasta2>) {
	chomp;
    if ($.>$start && $.< $end) {
        $seq .= $_;
    }
}

my $abase;
my $cbase;
my $gbase;
my $tbase;

$abase=($seq =~ tr/A//);
$cbase=($seq =~ tr/C//);
$gbase=($seq =~ tr/G//);
$tbase=($seq =~ tr/T//);

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

if ($cpg) {
	$temp_fasta = 'temp.fasta';
	open(TEMP, '>', $temp_fasta) or die "can't write to $temp_fasta: $!\n";
	print TEMP ">chr$chr\n";
	print TEMP "$seq\n";
	
	print "Running CpG Island analysis...\n";
	my $cpgicmd="perl CpGcluster.pl temp.fasta 50 1E-5 > CpGCluster.log";
	&forkExecWait($cpgicmd);
	print "Done\n";
	
	my $cpgi_file = "$wdir/temp.cpg";
	open my $cpgi, '<', $cpgi_file or die "can't open $file: $!";

	while (<$cpgi>) {
		push(@cpgi_index, $_);
	}
}

##############################################################################
#Count number of each nucleotide per bin (can be parallelized)
##############################################################################
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

##############################################################################
#Extract position column from the summary file
##############################################################################
print "Getting positions from summary file...\n";
my @POS;
#readline($index); #<-throws out summary header if it exists
while (<$index>) {
	push (@POS, (split(/\t/, $_))[1]);
}
print "Done\n";

##############################################################################
#Compare summary to reference and output pos/local subsequence/bin
#Update here to include adjacent sites
##############################################################################
print "Creating data file...\n";

if ($cpg) {
	
	print OUT "POS\tPAIR\tCPGI\n";

	foreach my $row (@POS) {
		my $hit=0;	
		foreach my $cpgi_int (@cpgi_index) {
			my @pair = split(/,/, $cpgi_int);
			my $start=$pair[0];
			my $end=$pair[1];
			
			if (($row >= $start) && ($row <= $end)) {
				$hit=1;
				last;
			}
		}

		print OUT "$row\t";
	    print OUT substr($seq, $row-$adj-1, $subseq);
		print OUT "\t";
		#print OUT ceil($row/$binwidth);
		#print OUT "\t";
		print OUT "$hit";
		print OUT "\n";
	}
} else {

	print OUT "POS\tPAIR\n";

	foreach my $row (@POS) {
		print OUT "$row\t";
	    print OUT substr($seq, $row-$adj-1, $subseq);
		#print OUT "\t";
		#print OUT ceil($row/$binwidth);
		print OUT "\n";
	}
}

print "Done\n";

##############################################################################
#Run selected R script
##############################################################################
print "Running R script...\n";
&forkExecWait($cmd);
print "Done. See images folder for output.\n";

##############################################################################
#Clean up temp files
##############################################################################
if ($cpg) {
	unlink $temp_fasta;
}
unlink $cpg_out;
unlink $bin_out;

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
# Hyun's subroutine--executes system commands
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

__END__
=head1 NAME

ref5.pl - Customized mysqldump utility

=head1 SYNOPSIS

        ref5.pl [OPTIONS]
        Options:
		--help			program documentation
		--chr			chromosome
		--mac			minor allele count
		--script		R script
		--b			binwidth
		--f			flank

=head1 OPTIONS

=over 8

=item B<--help>

Print a brief help message and exits.

=item B<--chr>

specify chromosome for analysis

=item B<--mac>

specify minor allele count of sites in existing summary file

=item B<--script>

specify which R script to run: 1 for counts, 2 for proportions

=item B<--b>

specify bin width for histograms (default is 100,000)

=item B<--f>

specify number of nucleotides in either direction from the variant to include in analysis
default includes only the next 5' nucleotide for CpG distinction

=back

=head1 DESCRIPTION

B<my-prog.pl> is doing something.

=head1 AUTHOR

=over

Jedidiah Carlson E<10> Department of Biostatistics E<10> University of Michigan

=cut