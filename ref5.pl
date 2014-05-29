#!/usr/local/bin/perl
use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;

##############################################################################
# SMAUG: Singleton Mutation Analysis Utility with Graphics
#
# Compares summary file from bcftools with reference genome and 
# outputs adjacent sites, then passes output to R script to plot distribution
#
# To-Do:
# -update directories to be consistent with pipe.pl
# -Parallel handling of fasta start/endpoints
# -Integrate with primary Perl pipeline
#
# Associated Files:
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
my $help=0;
my $man=0;
my $chr;
my $mac;
my $script;
my $binwidth=100000;
my $adj=0;

my $subseq=2;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

GetOptions ('chr=i'=> \$chr,
'mac=i'=> \$mac,
'script=i' => \$script,
'b=i' => \$binwidth,
'f=i' => \$adj,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;

if (!$chr | !$mac | !$script) {
	pod2usage("$0: Missing mandatory argument.");
}
##############################################################################
#Process mandatory inputs (need to change to elsif to prevent bad entries)
##############################################################################
my $dir;
if ($mac==1) {
	$dir = "singletons";
}
if ($mac==2) {
	$dir = "doubletons";
}

my $cmd;
my $args="$chr $dir $binwidth";
if ($script==1) {
	$cmd="Rscript count.R $args";
}
if ($script==2) {
	$cmd="Rscript prop.R $args";
}

my $sindex=$chr-1;
my $eindex=$chr;

##############################################################################
#Read in reference fasta file and summary file
##############################################################################
my $file = "/net/bipolar/jedidiah/human_g1k_v37.fasta";
open my $input, '<', $file or die "can't open $file: $!";

my $summ = "/net/bipolar/jedidiah/bcftools/summaries/$dir/all/chr$chr.$dir.summary.txt";
open my $index, '<', $summ or die "can't open $summ: $!";

##############################################################################
#initialize outputs
##############################################################################
my $cpg_out = 'cpg_out.txt';
open(OUT, '>', $cpg_out) or die "can't write to $cpg_out: $!\n";

my $bin_out = 'bin_out.txt';
open(BIN, '>', $bin_out) or die "can't write to $bin_out: $!\n";

my $fasta_out = 'temp.fasta';
open(TEMP, '>', $fasta_out) or die "can't write to $fasta_out: $!\n";

##############################################################################
#Find all headings from fasta file to define endpoints
##############################################################################
my @endpts;
while (<$input>) {
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
##############################################################################
#Extract only the position column from the summary file
##############################################################################
my @POS;
#readline($index); #<-throws out summary header if it exists
while (<$index>) {
	push (@POS, (split(/\t/, $_))[1]);
}

##############################################################################
#Re-read in reference fasta file and subset for selected chromosome
##############################################################################
my $file2 = 'human_g1k_v37.fasta';
open my $input2, '<', $file2 or die "can't open $file: $!";

my $seq;
while (<$input2>) {
	chomp;
    if ($.>$start && $.< $end) {
        $seq .= $_;
    }
}

print TEMP "$seq\n";

##############################################################################
#call CpGI script
##############################################################################
#my $cpgicmd="perl cpgi130.pl temp.fasta";
#&forkExecWait($cpgicmd);

##############################################################################
#Count number of each nucleotide per bin (can be parallelized)
##############################################################################
my $length=length($seq);
my $numbins=ceil($length/$binwidth);
my @A;
my @C;
my @G;
my @T;

for my $i (0 .. $numbins-1) {
	$A[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/A//);
    $C[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/C//);
    $G[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/G//);
    $T[$i]= (substr($seq, $i*$binwidth, $binwidth) =~ tr/T//);

	my $sum_at=$A[$i]+$T[$i];
	my $sum_cg=$C[$i]+$G[$i];

	print BIN "$sum_at\t$sum_cg\n";
}

#Verify output
#print "$length\n$numbins\n$A\n";
#foreach my $count (@A) {
#	print "$count\n";
#}

##############################################################################
#Compare summary to reference and output pos/local subsequence/bin
#Update here to include adjacent sites
##############################################################################
foreach my $row (@POS) {
    print OUT "$row\t";
    print OUT substr($seq, $row-$adj-1, $subseq);
	print OUT "\t";
	print OUT ceil($row/$binwidth);
    print OUT "\n";
}
##############################################################################
#Run selected R script
##############################################################################
#&forkExecWait($cmd);

##############################################################################
#Remove temp files
##############################################################################
#unlink $fasta_out;
#unlink $cpg_out;
#unlink $bin_out;

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