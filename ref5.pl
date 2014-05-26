#!/usr/local/bin/perl
use strict;
use warnings;
use POSIX;

###################################################
#Compares summary file from bcftools with 
#reference genome and outputs adjacent sites
#then passes output to R script to plot distribution
#
#To-Do:
#-code inputs as script arguments
#-Expand to more complex adjacent-site combinations
#-Parallel handling of fasta start/endpoints
#-Integrate with primary Perl pipeline
#test checkout
###################################################


#Prompt inputs
print "enter chromosome number (23 for X; 24 for Y)\n";

my $chr = <>;
chomp $chr;
print "1 for Singletons; 2 for Doubletons\n";

my $mac = <>;

print "Select R script:\n
1. Plot counts\n
2. Plot proportions\n";

my $script = <>;

#Process inputs (need to change to elsif to prevent bad entries)
my $dir;
if ($mac==1) {
	$dir = "singletons";
}
if ($mac==2) {
	$dir = "doubletons";
}

my $cmd;
my $args="$chr $dir";
if ($script==1) {
	$cmd="Rscript test.R $args";
}
if ($script==2) {
	$cmd="Rscript bin3.R $args";
}

my $sindex=$chr-1;
my $eindex=$chr;

#Read in reference fasta file and summary file
my $file = 'human_g1k_v37.fasta';
open my $input, '<', $file or die "can't open $file: $!";

my $summ = "/net/bipolar/jedidiah/bcftools/summaries/$dir/all/chr$chr.$dir.summary.txt";
open my $index, '<', $summ or die "can't open $summ: $!";

#initialize outputs
my $cpg_out = 'cpg_out.txt';
open(OUT, '>', $cpg_out) or die "can't write to $cpg_out: $!\n";

my $bin_out = 'bin_out.txt';
open(BIN, '>', $bin_out) or die "can't write to $bin_out: $!\n";

#Find all headings from fasta file to define endpoints
my @endpts;
while (<$input>) {
	if ($_=~ />/) {
		push (@endpts, $.);
	}
}

#Specify endpoints for selected chromosome (can be parallelized)
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

#Extract only the position column from the summary file
my @POS;
#readline($index); #<-throws out summary header if it exists
while (<$index>) {
	push (@POS, (split(/\t/, $_))[1]);
}

#Initialize temporary fasta out (passed to cpgi script)
my $fasta_out = 'temp.fasta';
open(TEMP, '>', $fasta_out) or die "can't write to $fasta_out: $!\n";

#Re-read in reference fasta file
my $file = 'human_g1k_v37.fasta';
open my $input2, '<', $file or die "can't open $file: $!";

#Subset fasta for selected chromosome and concatenate sequence to single string
my $seq;
while (<$input2>) {
	chomp;
    if ($.>$start && $.< $end) {
        $seq .= $_;
    }
}

print TEMP "$seq\n";

#call CpGI script
my $cpgicmd="perl cpgi130.pl temp.fasta";
&forkExecWait($cpgicmd);

#Initialize bin variables
my $length=length($seq);
my $binwidth=1000000;
my $numbins=ceil($length/$binwidth);
my @A;
my @C;
my @G;
my @T;

#Count number of each nucleotide per bin (can be parallelized)
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


#Compare summary to reference and output position, adjacent site, and bin number
foreach my $row (@POS) {
    print OUT "$row\t";
    print OUT substr($seq, $row-1, 2);
	print OUT "\t";
	print OUT ceil($row/$binwidth);
    print OUT "\n";
}

#Run selected R script
&forkExecWait($cmd);

#############################################
#Hyun's subroutine--executes system commands
#############################################
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
