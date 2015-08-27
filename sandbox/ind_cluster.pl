#!/usr/local/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Cwd;

my $subsetcmd="bcftools query -i 'N_ALT=1' -r 20 -o /net/bipolar/jedidiah/testpipe/vcfs/tmp.summary -O z /net/bipolar/jedidiah/testpipe/vcfs/merged.new.vcf.gz";
&forkExecWait($subsetcmd);

my $f_summ="net/bipolar/jedidiah/testpipe/vcfs/tmp.summary";
open my $summ, '<', $f_summ or die "can't open $f_summ: $!";

my $phet1=-1000;
my $nhet1=-1000;

my $phet2=-1000;
my $nhet2=-1000;

my $phet3=-1000;
my $nhet3=-1000;

my $phet4=-1000;
my $nhet4=-1000;

my $phet5=-1000;
my $nhet5=-1000;

my @prevs=(-1000, -1000, -1000, -1000, -1000);
my @nexts=(-1000, -1000, -1000, -1000, -1000);

my @scount=();
while(<$summ>){
	chomp;
	my @line=split(/\t/, $_);
	
	my $pos=$line[1];
	my $mac=$line[2];
	
	if($mac>=377){
		$prevs[4]=$nexts[4];
		$nexts[4]=$pos;
	}
	
	if($mac==1){
		if(($pos > $prevs[4]-1000 && $pos < $prevs[4]+1000) || ($pos > $nexts[4]-1000 && $pos < $nexts[4]+1000)){
			$scount++;
		}
	}


}

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