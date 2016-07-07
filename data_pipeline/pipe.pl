#!/usr/local/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Cwd;

################################################################################
# Singleton Analysis Pipeline
# 	-create initial summary files to be passed to ref5.pl
################################################################################

################################################################################
# initialize directories
# -$dir is location of program and input files
# -$projdir is location of output; assumed to contain /vcfs subfolder containing
#	at least 1 vcf file
#
# check if input files exist and process ped file if needed
################################################################################
my $dir=getcwd;

my $projdir="/net/bipolar/jedidiah/testpipe";
my $vcfloc="$projdir/vcfs";
make_path("$projdir/summaries/doubletons",
	"$projdir/summaries/cases",
	"$projdir/summaries/controls");

# my $summloc="$projdir/summaries";
my $summloc="$projdir/summaries/common";
my $summloc2="$projdir/summaries/doubletons";

my $bcftools="/net/bipolar/jedidiah/bcftools/bcftools";
my $vcftools="/net/bipolar/jedidiah/vcftools_0.1.10/bin/vcftools";

my $vcfin=1;
# my @vcfs = </net/bipolar/lockeae/freeze4/vcfs/anno/final/*.vcf.gz>;
my @vcfs;
my @chrindex=(1..22);
foreach my $chr (@chrindex){
	push(@vcfs,
		"/net/bipolar/lockeae/final_freeze/snps/vcfs/chr$chr/chr$chr.filtered.sites.vcf.gz");
}
# my @vcfs = </net/bipolar/lockeae/final_freeze/snps/vcfs/*.vcf.gz>;

if ($vcfin!=1) {
	print "Copying VCFs to project directory and updating headers...\n";
	foreach my $file (@vcfs) {
		my $filename=fileparse($file);
		# my $subfile = substr($filename, index($filename, 'chr'), index($filename, 'anno'));
		my $subfile = substr($filename,
			index($filename, 'chr'),
			index($filename, 'modified'));
		my $chr = substr($subfile, 0, index($subfile, '.'));
		#print "$chr\n";
		# my $tabix ="tabix -r newheader.txt $file > $vcfloc/$chr.anno.vcf.gz";
		# &forkExecWait($tabix);

		# my $tabix="tabix -p vcf $file"
		my $cpvcf="cp $file $vcfloc/$chr.vcf.gz";
		&forkExecWait($cpvcf);
	}
	print "Done\n";
}

my $script = 1;

################################################################################
#Produce summary files of singleton sites for cases, controls, and combined vcfs
################################################################################
if ($script==1){
	foreach my $chr (1..22) {
        # my $filename=fileparse($file);
        # my $path=dirname($file);
		# my $chr = substr($filename, 0, index($filename, '.'));
		# my $file = "/net/bipolar/jedidiah/testpipe/vcfs/merged.ma.vcf.gz"; # singletons, including multiallelic sites
		# my $cmd="bcftools query -i 'FILTER=\"PASS\"' -r $chr -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AB\t%INFO/AN\n' $file > $summloc/chr$chr.summary";

		my $file = "/net/bipolar/jedidiah/testpipe/vcfs/merged.new.vcf.gz"; # common variants
		my $cmd="bcftools query -i 'AC>10 && FILTER=\"PASS\"' -r $chr -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AB\t%INFO/AN\n' $file > $summloc/chr$chr.summary";
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
