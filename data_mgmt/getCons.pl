#!/usr/local/bin/perl

##############################################################################
# Script downloads and unzips phastCons46way conservation score data
##############################################################################

my $outfile = "/net/bipolar/jedidiah/mutation/reference_data/phastconsfiles.txt";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";


# for my $i (1 .. 22) {
	# print OUT "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/phastCons46way/primates/chr$i.phastCons46way.primates.wigFix.gz\n";
# }

# my $getcmd="wget -P /net/bipolar/jedidiah/mutation/reference_data/ -i $outfile";
# &forkExecWait($getcmd);

for my $i (1 .. 22) {
	my $unzipcmd="gunzip /net/bipolar/jedidiah/mutation/reference_data/chr$i.phastCons46way.primates.wigFix.gz";
	&forkExecWait($unzipcmd);
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