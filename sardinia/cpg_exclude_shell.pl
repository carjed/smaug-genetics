#!/usr/local/bin/perl

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Cwd;
use Benchmark;

my $parentdir="/net/bipolar/jedidiah";

my @covs=(3,5,8,10,12,15,20,30,40,80);

foreach(@covs){
	my $cpgcmd="perl cpg_exclude.pl --cov $_ &";
	&forkExecWait($cpgcmd);
}


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