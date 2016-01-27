#!/usr/local/bin/perl

##################################################################
# Utility script for calculating mean coverage from .bam files
# Usage:
# samtools mpileup 
#	/net/bipolar/jedidiah/sardinia/40x-in/28009_chr20_40x.bam | 
#	awk '{print $4}' | perl mean_coverage.pl
##################################################################

($num,$den)=(0,0);
while ($cov=<STDIN>) {
    $num=$num+$cov;
    $den++;
}
$cov=$num/$den;
print "Mean Coverage = $cov\n";