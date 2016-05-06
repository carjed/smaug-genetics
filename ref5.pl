#!/usr/local/bin/perl

##############################################################################
# SMAUG: Singleton Mutation Analysis Utility with Graphics
##############################################################################
# SMAUG uses extremely rare variants to visualize changes in mutation rates
# across the genome. The ref5.pl script takes bcftools summary files and
#
# Jedidiah Carlson
# Department of Bioinformatics
# The University of Michigan
##############################################################################

##############################################################################
#Initialize inputs, options, and defaults and define errors
##############################################################################
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

my $wdir=getcwd;
my $parentdir="/net/bipolar/jedidiah/mutation";

my $help=0;
my $man=0;
my $chr;
my $mac;
my $binwidth=100000;
my $adj=0;
my $data="full";
my $cpg='';
my $hot='';
my $mask_flag='';
my $rscript=0;
my @annoin;

GetOptions ('chr=s'=> \$chr,
'mac=i'=> \$mac,
'b=i' => \$binwidth,
'adj=i' => \$adj,
'data=s' => \$data,
'cpg' => \$cpg,
'hot' => \$hot,
'mf' => \$mask_flag,
'rs=i' => \$rscript,
'anno=s' => \@annoin,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

pod2usage(0) if $help;
pod2usage(-verbose => 2) if $man;

if (!$chr | !$mac) {
	pod2usage("$0: Missing mandatory argument.");
}

print "Local subsequence and CpG command entered simultaneously--overriding CpG analysis\n" if ($cpg && $adj!=0);

##############################################################################
#Process inputs
##############################################################################

make_path("$parentdir/images/chr$chr");
make_path("$parentdir/output/${subseq}bp_${bw}k_${macl}")
# make_path("$parentdir/output/${subseq}bp_${bw}k_${macl}")
my $imgdir="$parentdir/images/chr$chr";

@annoin=split(',',join(',',@annoin));
my @annos=qw(Intergenic Intron Nonsynonymous Exon Synonymous Utr3 Utr5 Upstream Downstream Stop_Gain Stop_Loss Start_Loss Essential_Splice_Site Normal_Splice_Site);
my @exonic=qw(Nonsynonymous Synonymous Stop_Gain Stop_Loss Start_Loss Essential_Splice_Site Normal_Splice_Site);
my @useannos;

if (@annoin) {
	print "Verifying selected annotations...\n";
	for (@annoin) {
		if ($_ ~~ @annos) {
			print "$_: valid\n";
			push(@useannos, $_);
		} else {
			print "$_: invalid annotation. See help file.\n";
		}
	}

	if (@useannos) {
		my $annodir=join('_',@useannos);
		make_path("$parentdir/images/chr$chr/anno/$annodir");
		$imgdir="$parentdir/images/chr$chr/anno/$annodir";
	}
}

print "Plots will be created in: $imgdir\n";

my $macl;
if ($mac==1) {
	$macl = "singletons";
}
if ($mac==2) {
	$macl = "doubletons";
}
if ($mac==11){
	$macl="common";
}

my $cpg_flag;
if ($cpg && $adj==0) {
	$cpg_flag="on";
} else {
	$cpg_flag="off";
}

my $hot_flag;
if ($hot) {
	$hot_flag="on";
} else {
	$hot_flag="off";
}

my $subseq=2;
if ($adj!=0) {
	$subseq = $adj*2+1;
}

my $bw=$binwidth/1000;

my $nextchr='';
if ($chr lt '22') {
	$nextchr=$chr+1;
} elsif ($chr eq '22') {
	$nextchr="X";
} else {
	$nextchr="Y";
}

##############################################################################
# Read in files and initialize outputs
# download hg37 from nih.gov if missing
##############################################################################
####OLD####my $f_summ = "/net/bipolar/jedidiah/bcftools/summaries/$macl/all/chr$chr.$macl.summary.txt";

my $f_summ = "/net/bipolar/jedidiah/testpipe/summaries/${macl}_${data}/chr$chr.summary"; #main line for full processing
# my $f_summ = "$parentdir/smaug-sandbox/scripts/human_chimp_chr10.summary";
# my $f_summ = "/net/bipolar/jedidiah/testpipe/summaries/chr10.common.summary";
# my $f_summ = "/net/bipolar/jedidiah/cur_vs_anc.summary";
# my $f_summ = "/net/bipolar/jedidiah/testpipe/vcfs/mask/summary/chr20.summary";
# my $f_summ = "/net/bipolar/jedidiah/testpipe/vcfs/chr20.rq.summary";
open my $summ, '<', $f_summ or die "can't open $f_summ: $!";

my $outfile = "$parentdir/output/${subseq}bp_${bw}k_${macl}/chr$chr.expanded.summary";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

if ($mac==1) {
	print OUT "CHR\tPOS\tREF\tALT\tDP\tAN\tANNO\t";
	# print OUT "CHR\tPOS\tREF\tALT\t";
} elsif ($mac==2) {
	print OUT "CHR\tPOS\tREF\tALT\tANNO\t";
}

my $bin_out = "$parentdir/output/${subseq}bp_${bw}k_${macl}/chr$chr.bin_out.txt";
open(BIN, '>', $bin_out) or die "can't write to $bin_out: $!\n";

my $seq=&getRef();

sub getRef{
	my $f_fasta;
	if($mask_flag){
		$f_fasta = "$parentdir/reference_data/human_g1k_v37.mask.fasta";
	} else {
		$f_fasta = "$parentdir/reference_data/human_g1k_v37.fasta";
	}

	if (-e $f_fasta) {
		print "Using reference genome: $f_fasta\n";
	} else {
		print "Reference genome not found in directory. Would you like to download one? (y/n): ";
		my $choice = <>;
		chomp $choice;
		if ($choice eq "y") {
			my $dlcmd="wget -P $parentdir/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz";
			&forkExecWait($dlcmd);
			my $unzipcmd="gunzip $parentdir/human_g1k_v37.fasta";
			&forkExecWait($unzipcmd);
		} else {
			die "Please upload an appropriate reference genome to $parentdir/reference_data/ \n";
		}
	}

	open my $fasta, '<', $f_fasta or die "can't open $f_fasta: $!";
	print "Getting reference sequence for chromosome $chr...\n";

	my $seq;
	while (<$fasta>) {
		chomp;
		if (/>$chr /../>$nextchr /) {
			next if />$chr / || />$nextchr /;
			$seq .=$_;
		}
	}

	return $seq;
}

my $altseq=$seq;
$altseq =~ tr/ACGT/TGCA/;

my $seqlength=length($seq);
print "seqlength: $seqlength\n";
print "Done\n";

##############################################################################
# Create index files for CpG islands if selected
##############################################################################

my @cpgi_index;
my $temp_fasta;

if ($cpg && $adj==0) {
	$temp_fasta = 'temp.fasta';
	open(TEMP, '>', $temp_fasta) or die "can't write to $temp_fasta: $!\n";
	print TEMP ">chr$chr\n";
	print TEMP "$seq\n";

	print "Running CpG Island analysis...\n";
	my $cpgicmd="perl CpGcluster.pl temp.fasta 50 1E-5 > CpGCluster.log";
	&forkExecWait($cpgicmd);
	print "Done\n";

	my $f_cpgi = "$wdir/temp.cpg";
	open my $cpgi, '<', $f_cpgi or die "can't open $f_cpgi: $!";

	while (<$cpgi>) {
		push(@cpgi_index, $_);
	}
}

##############################################################################
# Counts possible mutable sites per bin for 6 main categories
# and for local sequence analysis if selected
# -also returns GC content per bin
##############################################################################

&binCounts();

sub binCounts{
	print "Getting bin counts...\n";
	my $length=length($seq);
	print "seqlength: $length\n";
	my $numbins=ceil($length/$binwidth);
	my $bin;

	if ($adj>=1) {

		my @a= glob "{A,C,G,T}"x $subseq;
		my @b = (0) x (scalar @a);

		# &countSubSeq(@a, @b);

		print "Processing BIN1 header\n";
		print BIN "CHR\tAT\tCG\tprop_GC\tBIN\t";
		foreach my $tri (@a) {
			my $alttri=$tri;
			$alttri =~ tr/ACGT/TGCA/;
			# $alttri = reverse $alttri;
			my $min = minstr($tri, $alttri);
			my $max = reverse maxstr($tri, $alttri);
			# if ($tri !~ /^G|^T/) {
				# print BIN "$min($max)\t";
				print BIN "$tri\t";
			# }
		}
		print BIN "\n";

		print "Processing subseq counts per bin\n";
		for my $i (0 .. $numbins-1) {

			my $A = &nucCount($i, "A");
			my $C = &nucCount($i, "C");
			my $G = &nucCount($i, "G");
			my $T = &nucCount($i, "T");

			my @trinucs=(substr($seq, $i*$binwidth, $binwidth)=~/(?=([ACGT]{$subseq}))/g);
			#for(@trinucs){print "$_\n"};
			my %tri_count=();
			@tri_count{@a}=@b;
			$tri_count{$_}++ for @trinucs;

			my $sum_at=$A+$T;
			my $sum_cg=$C+$G;
			my $GC=0;
			if (($sum_at+$sum_cg)!=0) {
				$GC=($sum_cg/($sum_at+$sum_cg));
			}

			$bin=$i+1;

			print BIN "chr$chr\t$sum_at\t$sum_cg\t$GC\t$bin\t";
			#print BIN "$_:$tri_count{$_}\t" for sort keys(%tri_count);
			foreach my $count (sort keys %tri_count) {
				# if ($count !~ /N|^G|^T/) {
					my $altcount= $count;
					$altcount =~ tr/ACGT/TGCA/;
					$altcount = reverse $altcount;
					my $sum=$tri_count{$count}+$tri_count{$altcount};
					print BIN "$sum\t";
				# }
			}
			print BIN "\n";
		}
	} else {

		print BIN "CHR\tAT\tCG\tprop_GC\tBIN\n";

		for my $i (0 .. $numbins-1) {

			my $A = &nucCount($i, "A");
			my $C = &nucCount($i, "C");
			my $G = &nucCount($i, "G");
			my $T = &nucCount($i, "T");

			my $sum_at=$A+$T;
			my $sum_cg=$C+$G;
			my $GC=0;

			if (($sum_at+$sum_cg)!=0) {
				$GC=($sum_cg/($sum_at+$sum_cg));
			}

			$bin=$i+1;

			print BIN "chr$chr\t$sum_at\t$sum_cg\t$GC\t$bin\n";
		}
	}

	print "Done\n";
}

sub countSubSeq{
	my @seqs=shift;
	my @vec=shift;
	print "Counting subseqs in ref\n";
	my @trinucs=($seq=~/(?=([ACGT]{$subseq}))/g);
	my %tri_count=();
	@tri_count{@seqs}=@vec;
	$tri_count{$_}++ for @trinucs;
	print "Done\n";

	print "Processing BIN2\n";
	print BIN2 "SEQ\tCOUNT\n";
	foreach my $count (sort keys %tri_count) {
		if ($count !~ /N/) {
			print BIN2 "$count\t$tri_count{$count}\n";
		}
	}
	print "Done\n";
}

sub nucCount{
	my $binnum=shift;
	my $nuc=shift;

	my $out = () = substr($seq, $binnum*$binwidth, $binwidth) =~ /$nuc/g;
	return $out;
}

##############################################################################
# Output expanded summary file based on selected options
# -passed to R script along with bins file(s)
##############################################################################

print "Creating data file...\n";
my $start_time=new Benchmark;

my @POS;
my @NEWSUMM;
my @loci;
my $a_nu_start=0;
my $a_nu_start_cpg=0;
#readline($summ); #<-throws out summary header if it exists

if (@useannos) {
	print "subsetting sites by selected annotation(s)...\n";
	while (<$summ>) {
		next unless $_ =~ /^[^,]*$/;
		push (@POS, (split(/\t/, $_))[1]);
		foreach my $anno (@useannos) {
			if ($_ =~ /$anno/) {
				push (@NEWSUMM, $_);
			}
		}
	}
} else {
	print "No annotation subset selected. Using all data...\n";
	while (<$summ>) {
		next unless $_ =~ /^[^,]*$/;
		push (@POS, (split(/\t/, $_))[1]);
		push (@NEWSUMM, $_);
	}
}

if ($cpg && $adj==0) {

	print OUT "PAIR\tCPGI\tGC\n";

	foreach my $row (@NEWSUMM) {
		chomp $row;
		my @line=split(/\t/, $row);
		my $pos=$line[1];
		my $localseq = substr($seq, $pos-$adj-1, $subseq);
		my $cpgi = &getCpGI($pos);
		my $gcprop = &getGC($pos);

		print OUT "$row\t$localseq\t$cpgi\t$gcprop\n";
	}
} elsif ($hot) {

	&Hotspots;
	print OUT "SEQ\tALTSEQ\tGC\tDIST\n";
	#print OUT "CHR\tPOS\tREF\tALT\tDP\tAN\tANNO\tSEQ\tALTSEQ\tGC\tDIST\n";
	#print OUT "CHR\tPOS\tREF\tALT\tDP\tAN\tANNO\tSEQ\tALTSEQ\tGC\n";

	foreach my $row (@NEWSUMM) {
		chomp $row;
		my @line=split(/\t/, $row);
		my $pos=$line[1];
		my $localseq = substr($seq, $pos-$adj-1, $subseq);
		my $altlocalseq = reverse substr($altseq, $pos-$adj-1, $subseq);
		my $gcprop = &getGC($pos);
		my $distance = &getD2H($pos);

		print OUT "$row\t$localseq\t$altlocalseq\t$gcprop\t$distance\n";
		#print OUT "$row\t$localseq\t$altlocalseq\t$gcprop\n";
	}
} else {

	print OUT "SEQ\tALTSEQ\tGC\n";

	foreach my $row (@NEWSUMM) {
		chomp $row;
		my @line=split(/\t/, $row);
		my $pos=$line[1];
		my $localseq = substr($seq, $pos-$adj-1, $subseq);
		my $altlocalseq = reverse substr($altseq, $pos-$adj-1, $subseq);
		my $gcprop = &getGC($pos);

		if($localseq =~ /^[ACGT]+$/){
			print OUT "$row\t$localseq\t$altlocalseq\t$gcprop\n";
		}
	}
}

my $end_time=new Benchmark;
my $difference = timediff($end_time, $start_time);
print "Done. ";
print "Runtime: ", timestr($difference), "\n";

##############################################################################
#Run selected R script
##############################################################################

if($rscript==1){
	my $args="$chr $macl $binwidth $cpg_flag $outfile $adj $hot_flag $imgdir $bin_out";
	my $cmd="Rscript prop.R $args";

	print "Running R script...\n";
	&forkExecWait($cmd);
	print "Done. See images folder for output.\n";
}
##############################################################################
#Clean up temp files
##############################################################################

print "Cleaning up temp files...\n";

if ($cpg && $adj==0) {
	unlink $temp_fasta;
}

# unlink $outfile;
# unlink $bin_out;
# unlink $bin_out2;

my $plots_out="Rplots.pdf";
unlink $plots_out;

my $Rlog="R.log";
unlink $Rlog;

my $CpGlog="CpGCluster.log";
unlink $CpGlog;

my $CpGlog2="temp.cpg-log.txt";
unlink $CpGlog2;

my $CpGtemp="temp.cpg";
unlink $CpGtemp;

print "Done\n";

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

##############################################################################
# Recombination Hotspot Subroutine
# Returns 1 if site is inside or within 100bp of hotspot, 0 otherwise
##############################################################################
sub getD2H {
	my $site = shift;
	my $hit=0;

	my @first = split(/\t/, $loci[0]);
	my @last = split(/\t/, $loci[$#loci]);

	if ($site < $first[2]-100) {
		return $hit;
		last;
	} elsif ($site > $last[3]+100) {
		return $hit;
		last;
	}

	for my $i ($a_nu_start .. $#loci) {
		my $locus=$loci[$i];
		my @elements=split(/\t/, $locus);

		if ($site >=$elements[2]-100 && $site <= $elements[3]+100) {
			$hit=1;
			$a_nu_start=$i;
			return $hit;
			last;
		}
	}

	return $hit;
}

##############################################################################
# CpGI Subroutine
# Returns 1 if site is in CpG island, 0 otherwise
##############################################################################
sub getCpGI {
	my $site = shift;
	my $hit=0;

	my @first = split(/,/, $cpgi_index[0]);
	my @last = split(/,/, $cpgi_index[$#cpgi_index]);

	if ($site < $first[0]) {
		return $hit;
		last;
	} elsif ($site > $last[1]) {
		return $hit;
		last;
	}

	for my $i ($a_nu_start_cpg .. $#cpgi_index) {
		my $cpgi_int = $cpgi_index[$i];
		my @pair = split(/,/, $cpgi_int);

		if (($site >= $pair[0]) && ($site <= $pair[1])) {
			$hit=1;
			$a_nu_start_cpg=$i;
			return $hit;
			last;
		}
	}

	return $hit;
}

##############################################################################
# GC Content Subroutine
# Returns proportion of G or C nucleotides in reference genome for given
# flanking region of site
##############################################################################
sub getGC {
	my $site=shift;
	my $region_s=$site-($binwidth/2);
	#my $region_e=$site-$binwidth/2;
	my $region=substr($seq,$region_s,$binwidth);

	my $abase=($region =~ tr/A//);
	my $cbase=($region =~ tr/C//);
	my $gbase=($region =~ tr/G//);
	my $tbase=($region =~ tr/T//);

	my $gcsum=$cbase+$gbase;
	my $total=$abase+$cbase+$gbase+$tbase;
	my $gc_content=0.5;

	if ($total != 0) {
		$gc_content = $gcsum/$total;
	}
	return $gc_content;
}

# subgetExonic{
	# my$site
# }

##############################################################################
# Output expanded hotspot files
##############################################################################
sub Hotspots {
	my $f_hotspots = "$parentdir/reference_data/genetic_map/hotspots.txt"; #<-original hotspots input
	open my $hotspots, '<', $f_hotspots or die "can't open $f_hotspots: $!";

	my $hotspots_bed = "$parentdir/reference_data/genetic_map/hotspots.bed"; #<-initialize bed output
	open(HOTBED, '>', $hotspots_bed) or die "can't write to $hotspots_bed: $!\n";

	print "Analyzing recombination hotspots...\n";

	####### Read in original hotspots
	my $chrst;
	readline($hotspots); #<-throws out header
	while (<$hotspots>) {
		my @line=split(/\t/, $_);
		$chrst=$line[0];
		if ($chrst eq "chr$chr") {
			push (@loci, $_);
		}
	}

	####### Output hotspots in bed format
	foreach (@loci) {
		chomp;
		my @elements=split(/\t/, $_);
		my $chr=$elements[0];
		my $start=$elements[2];
		my $end=$elements[3];
		print HOTBED "$chr\t$start\t$end\n";
	}

	####### liftOver to hg19
	my $liftOvercmd="$parentdir/reference_data/liftOver $hotspots_bed $parentdir/reference_data/hg17ToHg19.over.chain $parentdir/reference_data/genetic_map/hotspots_hg19.bed unMapped";
	&forkExecWait($liftOvercmd);

	####### Read in liftOver output and per-site data
	my $f_new_hotspots = "$parentdir/reference_data/genetic_map/hotspots_hg19.bed";
	open my $new_hotspots, '<', $f_new_hotspots or die "can't open $f_new_hotspots: $!";

	my $f_site_data = "$parentdir/reference_data/genetic_map/genetic_map_GRCh37_chr$chr.txt";
	open my $site_data, '<', $f_site_data or die "can't open $f_site_data: $!";

	####### Initialize outputs
	my $hotspot_counts = 'hotspot_counts.txt';
	open(HOTCOUNT, '>', $hotspot_counts) or die "can't write to $hotspot_counts: $!\n";

	my $new_sites = 'new_sites.txt';
	open(NEWSITES, '>', $new_sites) or die "can't write to $new_sites: $!\n";

	my $bin_out3 = 'bin_out3.txt';
	open(BIN3, '>', $bin_out3) or die "can't write to $bin_out3: $!\n";

	my @newloci;
	while (<$new_hotspots>) {
		push (@newloci, $_);
	}

	print NEWSITES "Chromosome\tPOS\tRate\tMap\tStart\n";

	readline($site_data);
	while(<$site_data>) {
		chomp;
		my @data=split(/\t/, $_);
		my $site=$data[1];
		my $hotspot;

		foreach my $row (@newloci) {
			chomp;
			my @elements=split(/\t/, $row);

			my $start=$elements[1]-100;
			my $end=$elements[2]+100;
			my $length=$end-$start;
			my $center=$start+$length/2;

			if ($site >=$start && $site <= $end) {
				$hotspot=$start+100;
				last;
			} else {
				$hotspot=0;
			}
		}

		if ($hotspot!=0) {
			print NEWSITES "$_\t$hotspot\n";
		}
	}

	####### Process hg19 hotspots
	my $summind=0;
	my $oldend=0;
	my $rc_seq;
	print HOTCOUNT "Chromosome\tStart\tEnd\tCentre\tWidth(kb)\tD2E\tAT\tGC\tATCOUNT\tGCCOUNT\tAT_CG\tAT_GC\tAT_TA\tGC_AT\tGC_CG\tGC_TA\n";
	foreach (@newloci) {
		chomp;
		my @elements=split(/\t/, $_);

		my $start=$elements[1]-100;
		my $end=$elements[2]+100;
		my $length=$end-$start;
		my $center=$start+$length/2;

		my $d2e=$start-$oldend;
		$oldend=$end;

		#Count tri-nucleotide sequences within hotspots
		my $hotspot_seq=substr($seq,$start-1,$length);
		$rc_seq.="$hotspot_seq.";
		#print "$rc_seq\n";

		my $abase=($hotspot_seq =~ tr/A//);
		my $cbase=($hotspot_seq =~ tr/C//);
		my $gbase=($hotspot_seq =~ tr/G//);
		my $tbase=($hotspot_seq =~ tr/T//);

		my $gcsum=$cbase+$gbase;
		my $atsum=$abase+$tbase;

		print HOTCOUNT "$_\t$center\t$length\t$d2e\t$atsum\t$gcsum\t";

		my $atcount=0;
		my $gccount=0;

		my $AT_CG=0;
		my $AT_GC=0;
		my $AT_TA=0;

		my $GC_AT=0;
		my $GC_CG=0;
		my $GC_TA=0;

		for my $i ($summind .. $#NEWSUMM) {
			my $row = $NEWSUMM[$i];

			my @line=split(/\t/, $row);
			my $pos=$line[1];
			my $ref=$line[2];
			my $alt=$line[3];

			if ($pos >=$start && $pos <= $end) {
				#print "$pos:$ref\n";
				$atcount++ if ($ref eq "A" or $ref eq "T");
				$gccount++ if ($ref eq "C" or $ref eq "G");

				$AT_CG++ if ($ref eq "A" and $alt eq "C") or ($ref eq "T" and $alt eq "G");
				$AT_GC++ if ($ref eq "A" and $alt eq "G") or ($ref eq "T" and $alt eq "C");
				$AT_TA++ if ($ref eq "A" and $alt eq "T") or ($ref eq "T" and $alt eq "A");

				$GC_AT++ if ($ref eq "G" and $alt eq "A") or ($ref eq "C" and $alt eq "T");
				$GC_CG++ if ($ref eq "G" and $alt eq "C") or ($ref eq "C" and $alt eq "G");
				$GC_TA++ if ($ref eq "G" and $alt eq "T") or ($ref eq "C" and $alt eq "A");

				$summind=$i;
			}

			if ($pos > $end) {
				$summind=$i;
				last;
			}
		}

		print HOTCOUNT "$atcount\t$gccount\t$AT_CG\t$AT_GC\t$AT_TA\t$GC_AT\t$GC_CG\t$GC_TA\n";
	}

	my @a= glob "{A,C,G,T}"x $subseq;
	my @b = (0) x (scalar @a);

	my @trinucs=($rc_seq=~/(?=([ACGT]{$subseq}))/g);
	my %tri_count=();
	@tri_count{@a}=@b;
	$tri_count{$_}++ for @trinucs;


	print BIN3 "SEQ\tCOUNT\n";
	foreach my $count (sort keys %tri_count) {
		if ($count !~ /N/) {
			print BIN3 "$count\t$tri_count{$count}\n";
		}
	}
}

__END__
=head1 NAME

ref5.pl - SMAUG: Singleton Mutation Analysis Utility with Graphics

=head1 SYNOPSIS

        ref5.pl [OPTIONS]
        Options:
		--help			program documentation
		--chr			chromosome
		--mac			minor allele count
		--b			binwidth
		--adj			number of adjacent nucleotides
		--data			data subset to use
		--cpg			CpG site analysis?
		--hot			recombination hotspots?
		--anno			annotation(s)

=head1 OPTIONS

=over 8

=item B<--help>

Display this documentation

=item B<--chr>

MANDATORY: specify chromosome for analysis

=item B<--mac>

MANDATORY: specify minor allele count of sites in existing summary file

=item B<--b>

specify bin width for histograms (default is 100,000)

=item B<--adj>

specify number of adjacent nucleotides in either direction from the variant to include in analysis
default includes only the adjacent 3' nucleotide for CpG distinction

=item B<--data>

specify whether to use summaries from all singletons (full) or those that pass the strict filters (strict)

=item B<--cpg>

toggles extra analysis specific to CpG sites

=item B<--hot>

toggles extra analysis for distance to nearest recombination hotspot

=item B<--anno>

comma-separated list consisting of any of the following annotations:

Intergenic
Intron
Nonsynonymous
Exon
Synonymous
Utr3
Utr5
Upstream
Downstream
Stop_Gain
Stop_Loss
Start_Loss
Essential_Splice_Site
Normal_Splice_Site

=back

=head1 DESCRIPTION

B<ref5.pl> annotates summary files and counts motifs in genome over fixed-width windows

=head1 AUTHOR

=over

Jedidiah Carlson E<10> Department of Bioinformatics E<10> University of Michigan

=cut
