#!/usr/local/bin/perl

##############################################################################
# Script maps wigFix conservation scores to starch files
##############################################################################

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Cwd;

my $parentdir="/net/bipolar/jedidiah/mutation";
my $bedmap="/net/bipolar/jedidiah/bin/bedmap";

my $rootPath="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phastCons46way/primates/";
for my $i (1 .. 2){
	
	print "Mapping wigFix to starch for chromosome chr$i...\n";
	my $wigFn="chr$i.phastCons46way.primates.wigFix";
	print "$parentdir/reference_data/$wigFn\n";
	my $url="$rootPath/$wigFn.gz";
	# my $cmd="wget -qO- $url | gunzip - | /net/bipolar/jedidiah/bin/wig2starch - > $wigFn.starch";
	my $cmd="/net/bipolar/jedidiah/bin/wig2bed -i wig $parentdir/reference_data/$wigFn > $parentdir/reference_data/phastCons/$wigFn.bed";
	&forkExecWait($cmd);
	
	print "Mapping starch to 100kb...\n";
	my $mapcmd="$bedmap --echo --mean --chrom chr$i $parentdir/reference_data/genome.100kb.sorted.bed $parentdir/reference_data/phastCons/$wigFn.bed > $parentdir/reference_data/phastCons/phastCons46way.${i}.bed";
	&forkExecWait($mapcmd);
}

print "Done\n";

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