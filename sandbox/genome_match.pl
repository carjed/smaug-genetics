#!/usr/local/bin/perl

##############################################################################
# Script used to match positions where ancestral allele is same as current
# allele in reference sequences
##############################################################################

# use PDL; 
# use PDL::Char;                                                                  
# $PDL::SHARE=$PDL::SHARE; # keep stray warning quiet 

# my $source=PDL::Char->new("ATTCCGGG");                                          
# for my $str ( "ATTGCGGG", "ATACCGGC") {                                         
  # my $match =PDL::Char->new($str);                                              
  # my @diff=which($match!=$source)->list;                                        
  # print "@diff\n";                                                              
# }

my $chr=10;
my $nextchr=11;

my $a_fasta_file = "/net/1000g/hmkang/data/ref/human_ancestor_GRCh37_e59_hmk.fa";
open my $a_fasta, '<', $a_fasta_file or die "can't open $a_fasta_file: $!";

my $c_fasta_file = "/net/bipolar/jedidiah/human_g1k_v37.fasta";
open my $c_fasta, '<', $c_fasta_file or die "can't open $c_fasta_file: $!";

print "Getting reference sequence for chromosome $chr...\n";

my $anc_seq;
while (<$a_fasta>) {
	chomp;
	if (/>$chr/../>$nextchr/) {
		next if />$chr/ || />$nextchr/;
		$anc_seq .=$_;
	}
}
$anc_seq=uc($anc_seq);


my $cur_seq;
while (<$c_fasta>) {
	chomp;
	if (/>$chr/../>$nextchr/) {
		next if />$chr/ || />$nextchr/;
		$cur_seq .=$_;
	}
}
$cur_seq=uc($cur_seq);

# my $anc_altseq=$anc_seq;
# $anc_altseq =~ tr/ACGT/TGCA/;

# my $cur_seq="NNNACTACG";
# my $anc_seq="...ACTCcG";

my $outfile = "cur_vs_anc.summary";
open(OUT, '>', $outfile) or die "can't write to $outfile: $!\n";

for my $i (0 .. length($cur_seq)){
	my $c_char = substr($cur_seq, $i, 1);
	my $a_char =substr($anc_seq, $i, 1);
	# print "$c_char\t$a_char\n";
    if($c_char ne $a_char) {
		if($c_char =~ /[ACGT]/ && $a_char =~ /[ACGT]/){
			print OUT "$chr\t$i\t$a_char\t$c_char\n";
		}
	}
	
}	
	


# print "Done\n";