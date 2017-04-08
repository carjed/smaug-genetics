#!/usr/local/bin/perl

##############################################################################
# Run augment_summary.pl on all chromosomes--
# update local nucleotide (adj) and binwidth (b) parameters as necesary
#
# Can be modified to a slurm array job
##############################################################################

for my $i (1 .. 22) {
	my $cmd = "perl augment_summary.pl --chr $i &";
	# my $cmd = "perl augment_summary.pl --chr $i --mac 1 --adj 1 --b 1000000 --data full --motifs --binscheme band &";
	&forkExecWait($cmd);
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
