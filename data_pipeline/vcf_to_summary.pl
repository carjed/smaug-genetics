#!/usr/local/bin/perl

################################################################################
# Singleton Analysis Pipeline
# 	-create initial summary files to be passed to ref5.pl
################################################################################

use strict;
use warnings;
use POSIX;
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

# Default summary directory
my $outdir="$inputdir/summaries/${mac}";

# Copies per-chromosome VCFs from another directory to avoid overwriting
my $makecopy='';
if(defined($ARGV[0])){
  $makecopy = $ARGV[0];
}


# Specify project folder for VCFs
my $vcfdir="$inputdir/vcfs";
make_path($vcfdir);

################################################################################
# Copies original vcfs to project directory
################################################################################
if ($makecopy eq "copy") {
	my @rawvcfs = File::Find::Rule->file()
                            ->name("*.$rawvcfext")
                            ->maxdepth(3)
                            ->in($rawvcfdir);

	print "Processing VCFs with extension '$rawvcfext' from $rawvcfdir...\n";
	foreach my $rawvcf (@rawvcfs) {
		my $filename=fileparse($rawvcf);

    if($filename !~ /chrX/){
      my @parts = split(/\.vcf.gz/, $filename);
      my $basename = $parts[0];
      my $newvcf = "$vcfdir/$basename.ma.vcf.gz";
      # my $cpvcfcmd="cp $file $vcfdir/$filename";
      # print "Copying: $cpvcfcmd";
      # forkExecWait($cpvcfcmd);
      # print "Done\n";

      my $maparse="perl $parentdir/smaug-genetics/data_pipeline/ma_parse.pl --i $rawvcf --o $newvcf";
      print "Input file: $rawvcf\n";
      print "Writing to: $newvcf...";
      forkExecWait($maparse);
      print "Done\n";
    }
	}

	# my $concatcmd="perl $vcftoolsdir/perl/vcf-concat $vcfdir/*$rawvcfext | gzip -c > $vcfdir/merged.vcf.gz";
	# forkExecWait($concatcmd);

  print "Operation complete\n";
}

################################################################################
# Scans input directory for vcfs and outputs summary file to
# Get per-chromosome summary files from bcftools
################################################################################
my $script = 0;

if ($script==1){
  my @vcfs = File::Find::Rule->file()
                            ->name("*.vcf.gz")
                            ->maxdepth(1)
                            ->in($inputdir);

	foreach my $file (@vcfs) {

    if(!(-e "$file.tbi")){
      print "$file not indexed--indexing now...\n";
      my $tabixcmd = "tabix -p vcf $file";
      forkExecWait($tabixcmd);
      print "Done\n";
    }

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
    my $cmd = "$bcfquery  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AN\n' $file > $outdir/chr$chr.test.summary";
		forkExecWait($cmd);
	}
}
