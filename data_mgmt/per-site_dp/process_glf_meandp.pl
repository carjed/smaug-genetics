#!/usr/local/bin/perl

##############################################################################
# Scans summarized glf files (every 10bp; colnames: chr, pos, ref, dp) and
# outputs mean depth
##############################################################################

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Math::Round;
use Cwd;
use Benchmark;
use Tie::File;
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

my $chr;
my $ind;
my $help=0;
my $man=0;

# Set options and inputs
GetOptions ('chr=i'=> \$chr,
'ind=i' => \$ind,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;
# my $chr=1;

my $f_dirlist = "$parentdir/output/glf_depth/chr${chr}_glf_dirlist.txt";
open my $dirlist, '<', $f_dirlist or die "can't open $f_dirlist: $!";

my @dirs;
while(<$dirlist>){
  chomp;
  push @dirs, $_;
}

my $dir=$dirs[$ind-1];

my @path=split/\//, $dir;
my $chunk=$path[-1];

my @range = split/\./, $chunk;
my $start = roundup($range[0], 10);
my $end = roundup($range[1], 10)-10;

my %hash=();
my %hashn=();

for(my $i=$start; $i<=$end; $i+=10){
  $hashn{$i}=0;
}

my $outfile = "$parentdir/output/glf_depth/meandp/chr$chr.$chunk.txt";
print "$outfile\n";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

# glob ('/path/to/dir/*');
my @files = glob("$dir/*.dp");
my $i=0;
my $numfiles=scalar @files;
foreach my $file (@files) {
  # print $file . "\n";
  open my $sample, '<', $file or die "can't open $file: $!";

  while (<$sample>){
  	chomp;
  	my @line=split(/\t/, $_);
  	my $pos=$line[1];
  	# my $vals=join("\t", nearest(0.0001, @line[1 .. $#line]));
  	my $dp=$line[3];
  	$hash{$pos}+=$dp;
    $hashn{$pos}+=1;
  }

  # if($i%10==0){
  #   print "Finished $i of $numfiles samples\n";
  # }
  # $i++;
}

$hash{$_}=nearest(.01,$hash{$_}/$hashn{$_}) foreach (keys%hash);

print OUT "$_\t$hash{$_}\n" foreach (sort {$a <=> $b} (keys%hash));

print "Removing files in $dir/\n";
my $rmdpcmd="rm -f $dir/\*.dp";
&forkExecWait($rmdpcmd);
my $rmokcmd="rm -f $dir/samples.ok";
&forkExecWait($rmokcmd);

##############################################################################
# Subroutine rounds up to nearest roundto value
##############################################################################
sub roundup
{
  my $num = shift;
  my $roundto = shift || 1;

  return int(ceil($num/$roundto))*$roundto;
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
  } elsif ($kidpid==0) {
	  exec($cmd);
	  die "Cannot exec $cmd: $!";
  } else {
	  waitpid($kidpid,0);
  }
}
