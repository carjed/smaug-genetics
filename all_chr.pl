#!/usr/local/bin/perl

##############################################################################
# Run ref5.pl on all chromosomes--
# update local nucleotide (adj) and binwidth (b) parameters as necesary
# 
# Can be modified to a slurm array job
##############################################################################

for my $i (11 .. 22) {
	my $cmd = "perl ref5.pl --chr $i --mac 1 --adj 1 --b 100000 &";
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