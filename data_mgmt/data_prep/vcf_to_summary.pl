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
my $outdir="$parentdir/summaries";
make_path($outdir);

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
# Copies raw vcfs from original location to input directory
################################################################################
if ($makecopy eq "copy") {
	my @rawvcfs = File::Find::Rule->file()
                            ->name("*.$rawvcfext")
                            ->maxdepth(3)
                            ->in($rawvcfdir);

	# print STDERR "Processing VCFs with extension '$rawvcfext'\n from $rawvcfdir...\n";
	foreach my $rawvcf (@rawvcfs) {
		my $filename=fileparse($rawvcf);

    if($filename !~ /chrX/){
      my @parts = split(/\.vcf.gz/, $filename);
      my $basename = $parts[0];
      my @nameparts = split(/\./, $basename);
      my $i = $nameparts[0];
      $i =~ s/chr//g;

      my $newvcf = "$vcfdir/$basename.ma.aa.vcf.gz";
      my $ancestral = "$parentdir/reference_data/human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz";
      my $fasta = "$parentdir/reference_data//human_g1k_v37/chr$i.fasta.gz";

      # first command extracts singletons with filter PASS, including any that
      # occur in multiallelic sites
      my $maparse = "perl $relpath/ma_parse.pl --i $rawvcf";

      # second command fills ancestral allele to AA field
      my $aaparse = "perl $vcftoolsdir/perl/fill-aa -a $ancestral";

      # last command adds Motif and Category info fields
      my $infoparse = "perl $relpath/fill_motif.pl -a $fasta";

      # pipe commands and execute
      my $pipe = "$maparse | $aaparse | $infoparse | bgzip -c > $newvcf";
      print STDERR "Input file: $rawvcf\n";
      print STDERR "Writing to: $newvcf...\n";
      forkExecWait($pipe);
      print STDERR "Done\n";
    }
	}

  print STDERR "Operation complete\n";
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

  my $header = "\"CHR\\tPOS\\tREF\\tALT\\tAA\\tAN\\tMotif\\tCategory\"";
  my $headercmd = "echo $header > $outdir/full.summary";
  forkExecWait($headercmd);

	foreach my $file (@vcfs) {
    print STDERR "Getting summary for $file...\n";
    unless(-e "$file.tbi"){
      print STDERR "$file not indexed--indexing now...";
      my $tabixcmd = "tabix -p vcf $file";
      forkExecWait($tabixcmd);
      print STDERR "Done\n";
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
    my $outputcols = "'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\t%INFO/AN\t%INFO/Motif\t%INFO/Category\n'";
    # my $cmd = "$bcfquery -f $outputcols $file > $outdir/chr$chr.summary";
    my $cmd = "$bcfquery -f $outputcols $file >> $outdir/full.summary";
		forkExecWait($cmd);
	}
  print STDERR "Operation complete\n";
}
