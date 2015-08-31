#!/usr/local/bin/perl

##############################################################################
# Aggregate histone marks into 100kb windows
##############################################################################

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Cwd;


# my $bcftools="/net/bipolar/jedidiah/bcftools/bcftools";
# my $vcftools="/net/bipolar/jedidiah/vcftools_0.1.10/bin/vcftools";

my $refdir="/net/bipolar/jedidiah/mutation/reference_data";
my $binwidth="1000kb";
my $bedmap="/net/bipolar/jedidiah/bin/bedmap";

my @bedfiles=</net/bipolar/jedidiah/mutation/reference_data/histone_marks/*.bed>;
my @bedindex;
foreach my $file (@bedfiles){
	my $filename=basename($file, ".bed");
	my @strs=split(/\./, $filename);
	my $name= $strs[2];
	
	
	print "Binning $name ...\n";
	# print "$filename\n";
	
	my $cmd= "$bedmap --delim '\t' --echo --count --fraction-map 0.51 $refdir/genome.$binwidth.sorted.bed $file > $refdir/histone_marks/binned/$name.$binwidth.bed";
	
	&forkExecWait($cmd);
	
	print "Done\n";
	
}

my $cpgicmd="$bedmap --header --delim '\t' --echo --bases $refdir/genome.$binwidth.sorted.bed $refdir/cpg_islands_sorted.bed > $refdir/cpg_islands_$binwidth.bed";
&forkExecWait($cpgicmd);


my $exoncmd="$bedmap --delim '\t' --echo --bases $refdir/genome.$binwidth.sorted.bed $refdir/GRCh37_RefSeq_sorted.bed > /net/bipolar/jedidiah/mutation/reference_data/coding_exon_$binwidth.bed";
&forkExecWait($exoncmd);


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