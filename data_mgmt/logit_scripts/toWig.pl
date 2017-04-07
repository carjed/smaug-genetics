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
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));

my $config = LoadFile("$configpath/_config.yaml");

my $adj = $config->{adj};
my $mac = $config->{mac};
my $binw = $config->{binw};
my $data = $config->{data};
my $bin_scheme = $config->{bin_scheme};
my $parentdir = $config->{parentdir};
my $count_motifs = $config->{count_motifs};
my $expand_summ = $config->{expand_summ};

my $chr=$ARGV[0];

my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );

foreach my $categ (@categs) {
  my $outfile = "$parentdir/output/predicted/tracks/chr${chr}_$categ.wig";
  my $headercmd = "echo -e \"variableStep\\tchrom=chr${chr}\" > $outfile";
  &forkExecWait($headercmd);

  my $datacmd = "cut -f2,4 $parentdir/output/predicted/chr${chr}.$categ.txt >> $outfile";
  &forkExecWait($datacmd);
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
