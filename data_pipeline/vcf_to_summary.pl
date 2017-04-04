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
my $configpath = dirname($relpath);

my $config = LoadFile("$configpath/_config.yaml");

print "Script will run with the following parameters:\n";
for (sort keys %{$config}) {
    say "$_: $config->{$_}";
}

my $adj = $config->{adj};
my $mac = $config->{mac};
my $binw = $config->{binw};
my $data = $config->{data};
my $bin_scheme = $config->{bin_scheme};
my $parentdir = $config->{parentdir};
my $count_motifs = $config->{count_motifs};
my $expand_summ = $config->{expand_summ};
my $inputdir = $config->{inputdir};

################################################################################
# Singleton Analysis Pipeline
# 	-create initial summary files to be passed to ref5.pl
################################################################################
my $help=0;
my $man=0;

### Mandatory inputs
# Default vcf to load
my $invcf = "$inputdir/vcfs/merged.ma.vcf.gz"; # singletons, including multiallelic sites
# my $invcf = "/net/bipolar/jedidiah/testpipe/vcfs/merged.new.vcf.gz"; # common variants

# Default summary directory
my $outdir="$inputdir/summaries/${mac}_full";

### Optional inputs
# If --common option specified, outputs summary for common variants (MAC>10)
my $common;

# Copies per-chromosome VCFs from another directory to avoid overwriting
my $makecopy;

# Specify project folder for VCFs
my $vcfdir="$inputdir/vcfs";

GetOptions ('invcf=s'=> \$invcf,
'outdir=s'=> \$outdir,
'common' => \$common,
'makecopy' => \$makecopy,
'vcfdir=s' => \$vcfdir,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

################################################################################
# Copies original vcfs to project directory
# Must modify hard-coded paths for original vcfs
################################################################################
if ($makecopy) {
	my @vcfs;
	my @chrindex=(1..22);
	foreach my $chr (@chrindex){
		push(@vcfs,
			"/net/bipolar/lockeae/final_freeze/snps/vcfs/chr$chr/chr$chr.filtered.sites.vcf.gz");
	}

	print "Copying VCFs to project directory...\n";
	foreach my $file (@vcfs) {
		my $filename=fileparse($file);
		# my $subfile = substr($filename, index($filename, 'chr'), index($filename, 'anno'));
		my $subfile = substr($filename,
			index($filename, 'chr'),
			index($filename, 'modified'));
		my $chr = substr($subfile, 0, index($subfile, '.'));
		#print "$chr\n";
		# my $tabix ="tabix -r newheader.txt $file > $vcfdir/$chr.anno.vcf.gz";
		# &forkExecWait($tabix);

		# my $tabix="tabix -p vcf $file"
		my $cpvcf="cp $file $vcfdir/$chr.vcf.gz";
		&forkExecWait($cpvcf);
	}

	my $concatcmd="perl /net/bipolar/jedidiah/vcftools_0.1.10/perl/vcf-concat $vcfdir/chr*.vcf.gz | gzip -c > $vcfdir/merged.vcf.gz";
	&forkExecWait($concatcmd);

	my $maparse="perl ";

	print "Done\n";
}

################################################################################
# Get per-chromosome summary files from bcftools
################################################################################
my $script = 1;

if ($script==1){
	foreach my $chr (1..22) {
        # my $filename=fileparse($file);
        # my $path=dirname($file);
		# my $chr = substr($filename, 0, index($filename, '.'));
		my $cmd;
		if(common){
			$cmd="bcftools query -i 'AC>10 && FILTER=\"PASS\"' -r $chr -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AB\t%INFO/AN\n' $invcf > $outdir/chr$chr.common.summary";
		} else {
			$cmd="bcftools query -i 'FILTER=\"PASS\"' -r $chr -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AB\t%INFO/AN\n' $invcf > $outdir/chr$chr.summary";
		}

		&forkExecWait($cmd);
	}
}

################################################################################
# Fork-exec-wait
################################################################################
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
