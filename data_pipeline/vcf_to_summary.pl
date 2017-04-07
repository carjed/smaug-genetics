#!/usr/local/bin/perl
use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path qw(make_path);
use File::Find::Rule;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Cwd;
use Benchmark;
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname($relpath);
my $config = LoadFile("$configpath/_config.yaml");

my $mac = $config->{mac};
my $parentdir = $config->{parentdir};
my $inputdir = $config->{inputdir};
my $vcftoolsdir = $config->{vcftoolsdir};
my $rawvcfdir = $config->{rawvcfdir};
my $rawvcfext = $config->{rawvcfext};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef);

################################################################################
# Singleton Analysis Pipeline
# 	-create initial summary files to be passed to ref5.pl
################################################################################
my $help=0;
my $man=0;

# Default summary directory
my $outdir="$inputdir/summaries/${mac}";

# Copies per-chromosome VCFs from another directory to avoid overwriting
my $makecopy;

# Specify project folder for VCFs
my $vcfdir="$inputdir/vcfs";
make_path($vcfdir);

GetOptions ('invcf=s'=> \$invcf,
'outdir=s'=> \$outdir,
'makecopy' => \$makecopy,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

################################################################################
# Copies original vcfs to project directory
################################################################################
if ($makecopy) {
	my @rawvcfs = File::Find::Rule->file()
                            ->name("*.$rawvcfext")
                            ->maxdepth(3)
                            ->in($rawvcfdir);

	print "Copying VCFs to project directory...\n";
	foreach my $file (@rawvcfs) {
		my $filename=fileparse($file);

    if(!($filename !~ /chrX/)){
      my @parts = split(/\.vcf.gz/, $filename);
      my $basename = $parts[0];

      my $cpvcf="cp $file $vcfdir/$filename";
      print "sh: $cpvcf\n";
      forkExecWait($cpvcf);

      my $maparse="perl ./ma_parse.pl --i $vcfdir/$filename --o $vcfdir/$basename.ma.vcf.gz";
      forkExecWait($maparse);
    }
	}

	# my $concatcmd="perl $vcftoolsdir/perl/vcf-concat $vcfdir/*$rawvcfext | gzip -c > $vcfdir/merged.vcf.gz";
	# forkExecWait($concatcmd);


	print "Done\n";
}

################################################################################
# Scans input directory for vcfs and outputs summary file to
# Get per-chromosome summary files from bcftools
################################################################################
my $script = 1;

if ($script==1){
  my @vcfs = File::Find::Rule->file()
                            ->name("*.vcf.gz")
                            ->maxdepth(1)
                            ->in($inputdir);

	foreach my $file (@vcfs) {
    my $filename=fileparse($file);
    my @parts = split(/\./, $filename);
    my $chr = $parts[0];
    $chr =~ s/chr//;
    my $bcfquery;
		if($mac eq "common"){
			$bcfquery = "bcftools query -i 'AC>10 && FILTER=\"PASS\"' -r $chr";
		} elsif ($mac eq "singletons"){
			$bcfquery = "bcftools query -i 'AC=1 && FILTER=\"PASS\"' -r $chr";
		}
    my $cmd = "$bcfquery  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AN\n' $file > $outdir/chr$chr.summary";
		forkExecWait($cmd);
	}
}
