#!/usr/local/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Cwd;

##########################################################################################
# Singleton Analysis Pipeline
# 
# This series of scripts subsets the singletons from full vcfs, separates cases/controls,
# produces per-subject statistics, calculates the singleton TS/TV ratio by annotation,
# and counts the singletons per gene for cases and controls
#
# To-Do:
# -implement various R scripts for processing .ped file, plotting data, etc.
# -implement the mutation analysis script (ref5.pl + R scripts)
# -expand to doubletons/all variants
# -improve directory structure of project folder
# -fix TS/TV script
# -pass variables to R scripts
#
# Required files:
# -process_ped.R
# -singleton_vs_cov.R
# -.ped file; produces the following:
# 	-ped_out.ped; used for plotting script
# 	-list of cases, one per line (cases.txt)
# 	-list of controls, one per line (controls.txt)
# 	-full list of subjects, one per line (subjects.txt)
#
##########################################################################################

##########################################################################################
# initialize directories
# -$dir is location of program and input files
# -$projdir is location of output; assumed to contain /vcfs subfolder containing 
#	at least 1 vcf file
#
# check if input files exist and process ped file if needed
##########################################################################################
my $dir=getcwd;
my $casefile="$dir/cases.txt";
my $controlfile="$dir/controls.txt";
my $subjfile="$dir/subjects.txt";
my $ped="$dir/freeze4.20131216.v11.ped";

if (-e $ped) {
	print "using .ped file: $ped\n";
} else {
	die "ped file does not exist!\n";
}

if (-e $casefile && -e $controlfile && -e $subjfile) {
	print "using case file: $casefile\n";
	print "using control file: $controlfile\n";
	print "using subject file: $subjfile\n";
} else {
	print "Missing one of the following files:\n
		   $casefile\n
		   $controlfile\n
		   $subjfile\n
		   Creating files...\n";
	my $cmd="Rscript process_ped.R";
	&forkExecWait($cmd);
	print "Done\n";
}	

my $projdir="/net/bipolar/jedidiah/testpipe";
my $vcfloc="$projdir/vcfs";
make_path("$projdir/summaries/doubletons", "$projdir/summaries/cases", "$projdir/summaries/controls");
# my $summloc="$projdir/summaries";
my $summloc="$projdir/summaries";
my $summloc2="$projdir/summaries/doubletons";

my $bcftools="/net/bipolar/jedidiah/bcftools/bcftools";
my $vcftools="/net/bipolar/jedidiah/vcftools_0.1.10/bin/vcftools";

my $vcfin=1;
# my @vcfs = </net/bipolar/lockeae/freeze4/vcfs/anno/final/*.vcf.gz>;
my @vcfs;
my @chrindex=(1..22);
foreach my $chr (@chrindex){
	push(@vcfs, "/net/bipolar/lockeae/final_freeze/snps/vcfs/chr$chr/chr$chr.filtered.sites.modified.vcf.gz");
}
# my @vcfs = </net/bipolar/lockeae/final_freeze/snps/vcfs/*.vcf.gz>;

if ($vcfin!=1) {
	print "Copying VCFs to project directory and updating headers...\n";
	foreach my $file (@vcfs) {
		my $filename=fileparse($file);
		# my $subfile = substr($filename, index($filename, 'chr'), index($filename, 'anno'));
		my $subfile = substr($filename, index($filename, 'chr'), index($filename, 'modified'));
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

##########################################################################################
#Display menu and prompt input
##########################################################################################
print "Select script:\n
1. Create Singleton VCFs\n
2. Subset Cases/Controls\n
3. Obtain summary files\n
----------------------------\n
4. Per-subject singleton counts\n
5. TS/TV Summaries\n
6. Count singletons per gene\n
7. Processes per-subject singleton counts\n
8. Create Doubleton VCFs\n
9. Obtain doubleton summary files\n";
my $script = <>;

##########################################################################################
# Subset VCFs into singleton-only files
# cases and controls are matched for equal number per defined $subjfile;
# "-s" flag can be ommitted from bcftools script if we aren't looking at phenotype info
# (e.g. for mutation spectrum analysis)
# can also update process_ped.R to choose smart vs. simple subsetting
##########################################################################################
if ($script==1){
    my @files = <$vcfloc/*.vcf.gz>;
	make_path("$vcfloc/singletons");
	my $path="$vcfloc/singletons";
	
	print "verifying VCFs are indexed...\n";
	
    foreach my $file (@files) {
		#verify vcfs are indexed	
		my $tbi="$file.tbi";
		if (-e $tbi) {
			print "$file is already indexed\n";
		} else {
			my $cmd="tabix -p vcf $file";
			&forkExecWait($cmd);
			print "$file: done\n";
		}	
	
		my $filename=fileparse($file);
		my $chr = substr($filename, 0, index($filename, '.'));
		my $outfile = "$path/$chr.singletons.vcf.gz";
	
		my $cmd="bcftools query -i AC=1 -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AN\t.\n' $file > $summloc/$chr.summary &";
		&forkExecWait($cmd);
	
	
		# if (-e $outfile) {
			# print "$outfile already exists\n";
		# } else {
			# my $cmd="$bcftools view -x -c 1 -C 1 -s $subjfile -o $outfile -O z $file &";
			# my $cmd="bcftools query -i AC=1 -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' -o $outfile -O z $file &";
			# &forkExecWait($cmd);
		# }
	}
	
	print "Files are being processed in the background. Monitor progress via top command.\n";
}

##########################################################################################
# doubletons
##########################################################################################
if ($script==8){
    my @files = <$vcfloc/*.vcf.gz>;
	make_path("$vcfloc/doubletons");
	my $path="$vcfloc/doubletons";
	
	print "verifying VCFs are indexed...\n";
	
    foreach my $file (@files) {
		#verify vcfs are indexed	
		my $tbi="$file.tbi";
		if (-e $tbi) {
			print "$file is already indexed\n";
		} else {
			my $cmd="tabix -p vcf $file";
			&forkExecWait($cmd);
			print "$file: done\n";
		}	
	
		my $filename=fileparse($file);
		my $chr = substr($filename, 0, index($filename, '.'));
		my $outfile = "$path/$chr.doubletons.anno.vcf.gz";
		
		if (-e $outfile) {
			print "$outfile already exists\n";
		} else {
			my $cmd="$bcftools view -x -c 2 -C 2 -G -o $outfile -O z $file &";
			&forkExecWait($cmd);
		}
	}
	
	print "Files are being processed in the background. Monitor progress via top command.\n";
}

##########################################################################################
# Doubleton summaries
##########################################################################################
if ($script==9){
	my @files = <$vcfloc/doubletons/*.vcf.gz>;
	foreach my $file (@files) {
        my $filename=fileparse($file);
        my $path=dirname($file);
		my $chr = substr($filename, 0, index($filename, '.'));
		my $cmd="$bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ANNO\n' $file > $summloc2/$chr.summary &";
		&forkExecWait($cmd);
	}
}

##########################################################################################
# Subset Cases/Controls from singleton vcfs
#
# note "-G" flag in bcftools command--removes individual genotype info for smaller 
# file size and faster summaries, but prevents the use of case/control vcfs 
# for individual analysis
##########################################################################################
if ($script==2){
    my @files = <$vcfloc/singletons/*.vcf.gz>;
	make_path("$vcfloc/singletons/cases", "$vcfloc/singletons/controls");
	my $casepath = "$vcfloc/singletons/cases";
	my $controlpath = "$vcfloc/singletons/controls";
	
    foreach my $file (@files) {
		my $filename=fileparse($file);
		my $chr = substr($filename, 0, index($filename, '.'));
		my $outcase="$casepath/$chr.singletons.cases.vcf.gz";
		my $outcontrol="$controlpath/$chr.singletons.controls.vcf.gz";
		
		my $cmd1="$bcftools view -x -G -s $casefile -o $outcase -O z $file &";
		my $cmd2="$bcftools view -x -G -s $controlfile -o $outcontrol -O z $file &";
        &forkExecWait($cmd1);
		&forkExecWait($cmd2);
    }
}

##########################################################################################
#Produce summary files of singleton sites for cases, controls, and combined vcfs
##########################################################################################
if ($script==3){
	my @files = <$vcfloc/singletons/*.vcf.gz>;
	foreach my $file (@files) {
        my $filename=fileparse($file);
        my $path=dirname($file);
		my $chr = substr($filename, 0, index($filename, '.'));
		my $cmd="$bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AN\t%INFO/ANNO\n' $file > $summloc/$chr.summary &";
		&forkExecWait($cmd);
	}
	
	# my @cases = <$vcfloc/singletons/cases/*.vcf.gz>;
	# foreach my $file (@cases) {
        # my $filename=fileparse($file);
        # my $path=dirname($file);
		# my $chr = substr($filename, 0, index($filename, '.'));
		# my $cmd="$bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AN\t%INFO/ANNO\n' $file > $summloc/cases/$chr.cases.summary &";
		# &forkExecWait($cmd);
	# }
	
	# my @controls = <$vcfloc/singletons/controls/*.vcf.gz>;
	# foreach my $file (@controls) {
        # my $filename=fileparse($file);
        # my $path=dirname($file);
		# my $chr = substr($filename, 0, index($filename, '.'));
		# my $cmd="$bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AN\t%INFO/ANNO\n' $file > $summloc/controls/$chr.controls.summary &";
		# &forkExecWait($cmd);
	# }	
}

##########################################################################################
#Per-subject singleton counts
##########################################################################################
if ($script==4){
	my @files = </net/bipolar/lockeae/final_freeze/snps/vcfs/chr*/chr*.filtered.modified.vcf.gz>;
	make_path("$projdir/singletons");
	my $singdir="$projdir/singletons";
	
	foreach my $file (@files) {
		my $filename=fileparse($file);
		my $chr = substr($filename, 0, index($filename, '.'));
		my $cmd= "$vcftools --gzvcf $file --singletons --out $singdir/$chr";
		&forkExecWait($cmd);
	}
}


###################
#Processes per-subject singleton counts
###################
if ($script==7){
	make_path("$projdir/singletons");
	my $singdir="$projdir/singletons";
	
	my $cmd2="cat $singdir/*.singletons | cut -f5 | sort | uniq -c >> $singdir/full.singletons.count";
	&forkExecWait($cmd2);
	
	my $cmd3="Rscript singleton_vs_cov.R";
	&forkExecWait($cmd3);
}

##########################################################################################
#TS/TV (must already have summary files)
##########################################################################################
if ($script==5){
	my @annos = ("Intergenic", 
				 "Intron", 
				 "Nonsynonymous", 
				 "Exon", 
				 "Synonymous", 
				 "Utr3", 
				 "Utr5", 
				 "Upstream", 
				 "Downstream", 
				 "Stop_Gain", 
				 "Stop_Loss", 
				 "Start_Loss", 
				 "Essential_Splice_Site", 
				 "Normal_Splice_Site");
	my @files = <$summloc/*.summary>;
	
	foreach my $anno (@annos) {
		foreach my $file (@files) {
			my $filename=fileparse($file);
        	my $path=dirname($file);
			my $cmd="grep $anno $file | awk '($3 ~ /A/ && $4 ~ /C/) || ($3 ~ /A/ && $4 ~ /T/) || ($3 ~ /C/ && $4 ~ /A/) || ($3 ~ /C/ && $4 ~ /G/) || ($3 ~ /G/ && $4 ~ /C/) || ($3 ~ /G/ && $4 ~ /T/) || ($3 ~ /T/ && $4 ~ /A/) || ($3 ~ /T/ && $4 ~ /G/)' - > $file.$anno.tv.txt";
			my $cmd2="grep $anno $file | awk '($3 ~ /A/ && $4 ~ /G/) || ($3 ~ /G/ && $4 ~ /A/) || ($3 ~ /C/ && $4 ~ /T/) || ($3 ~ /T/ && $4 ~ /C/)'  - > $file.$anno.ts.txt";
			&forkExecWait($cmd);
		}
	}
}

##########################################################################################
#Count singletons per gene
#requires summary files
##########################################################################################
if ($script==6){
	my @annos = ("Intergenic", 
				 "Intron", 
				 "Nonsynonymous", 
				 "Exon", 
				 "Synonymous", 
				 "Utr3", 
				 "Utr5", 
				 "Upstream", 
				 "Downstream", 
				 "Stop_Gain", 
				 "Stop_Loss", 
				 "Start_Loss", 
				 "Essential_Splice_Site", 
				 "Normal_Splice_Site");
	my @files = <$summloc/*.summary>;

	print "Splitting files by annotation\n";
	foreach my $anno (@annos) {
		foreach my $file (@files) {
			my $filename=fileparse($file);
        	my $path=dirname($file);
			my $chr = substr($filename, 0, index($filename, '.'));
			my $cmd="grep $anno $file > $path/$anno.$chr.summary";
			&forkExecWait($cmd);
		}
	}
	print "Done\n";
	
	print "Combining exonic variants and splice sites\n";
	my $cmd2="cat $summloc/N* $summloc/St* $summloc/Es* > $summloc/exome.summary";
	&forkExecWait($cmd2);
	print "Done\n";
	
	print "Counting hits per gene\n";
	my $cmd3="cut -f5 $summloc/exome.summary | cut -d':' -f2 |sort | uniq -c >  $summloc/exome_genes.txt";
	&forkExecWait($cmd3);
	print "Done. See file at: \n";
}

##########################################################################################
#Hyun's subroutine--executes system commands
##########################################################################################
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
