#!/usr/local/bin/perl

##############################################################################
# Script concatenates chromosome vcfs from GotCloud output
##############################################################################

my @vcfs;

my $root="/net/bipolar/jedidiah/sardinia/80x-gw/beagle";

for my $i (1 .. 22) {
	push (@vcfs, "$root/chr$i/chr$i.filtered.PASS.beagled.vcf.gz");
}

my $input;
foreach (@vcfs) {
	#chomp;
	$input .= "$_ ";
}

#print "$input\n";

my $cmd = "/net/bipolar/jedidiah/bcftools/bcftools concat -O z $input > /net/bipolar/jedidiah/sardinia/80x_full.vcf.gz";

&forkExecWait($cmd);

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