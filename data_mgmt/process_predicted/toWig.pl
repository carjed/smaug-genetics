#!/usr/local/bin/perl

use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Path qw(make_path);
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

my $parentdir = $config->{parentdir};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

my $chr=$ARGV[0];

my @categs = qw( AT_CG AT_GC AT_TA GC_AT GC_CG GC_TA );

foreach my $categ (@categs) {
  my $outfile = "$parentdir/output/predicted/tracks/chr${chr}_$categ.wig";
  my $headercmd = "echo -e \"variableStep\\tchrom=chr${chr}\" > $outfile";
  forkExecWait($headercmd);

  my $datacmd = "cut -f2,4 $parentdir/output/predicted/chr${chr}.$categ.txt >> $outfile";
  forkExecWait($datacmd);
}
