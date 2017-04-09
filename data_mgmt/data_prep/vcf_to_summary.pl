#!/usr/local/bin/perl

################################################################################
# Get initial summary information from input vcfs
#
# To avoid the possibility of modifying raw vcf data, includes option to copy
# raw vcfs to a working project directory
################################################################################
use strict;
use warnings;
use POSIX;
use File::Basename;
use File::Path qw(make_path);
use File::Find::Rule;
use FindBin;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
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

# Running script with argument 'copy'
# copies VCFs from another directory to avoid overwriting raw data
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
      my @nameparts = split(/\./, $basename);
      my $i = $nameparts[0];
      $i =~ s/"chr"//g;
      my $newvcf = "$vcfdir/$basename.ma.aa.vcf.gz";

      my $maparse="perl ./ma_parse.pl --i $rawvcf | perl $vcftoolsdir/perl/fill-aa -a $parentdir/reference_data/human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz | bgzip -c > $newvcf";
      print "Input file: $rawvcf\n";
      print "Writing to: $newvcf...";
      forkExecWait($maparse);
      print "Done\n";
    }
	}

  print "Operation complete\n";
}

################################################################################
# Scans input directory for vcfs and outputs summary file to get per-chromosome
# summary files from bcftools. If running without the 'copy' argument, assumes
# vcfs in the input directory have already been parsed for multiallelic sites
################################################################################
my $script = 1;

if ($script==1){
  my @vcfs = File::Find::Rule->file()
                            ->name("*.vcf.gz")
                            ->maxdepth(1)
                            ->in($vcfdir);

	foreach my $file (@vcfs) {
    print "Getting summary for $file...\n";
    unless(-e "$file.tbi"){
      print "$file not indexed--indexing now...";
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
    my $cmd = "$bcfquery  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\t%INFO/AN\n' $file > $outdir/chr$chr.summary";
		forkExecWait($cmd);
	}
  print "Operation complete\n";
}
