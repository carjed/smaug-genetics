#!/usr/local/bin/perl

##############################################################################
# Script for filling ancestral allele in BRIDGES vcfs
##############################################################################

for my $i (1 .. 22) {
	# my $refcmd = "cat /net/bipolar/jedidiah/mutation/reference_data/human_ancestor_GRCh37_e59/human_ancestor_$i.fa  | sed 's,>.*,>$i,' | bgzip -c > /net/bipolar/jedidiah/mutation/reference_data/human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz";
	# &forkExecWait($refcmd);
	
	# my $indexcmd="samtools faidx /net/bipolar/jedidiah/mutation/reference_data/human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz";
	# &forkExecWait($indexcmd);
	
	# my $processcmd="zcat /net/bipolar/jedidiah/testpipe/vcfs/chr$i.vcf.gz | perl /net/bipolar/jedidiah/vcftools_0.1.10/perl/fill-aa -a /net/bipolar/jedidiah/mutation/reference_data/human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz | bgzip -c > /net/bipolar/jedidiah/testpipe/vcfs/ancestral/chr$i.aa.vcf.gz";
	my $processcmd="zcat /net/bipolar/lockeae/final_freeze/snps/vcfs/chr$i/chr$i.filtered.modified.vcf.gz | perl /net/bipolar/jedidiah/vcftools_0.1.10/perl/fill-aa -a /net/bipolar/jedidiah/mutation/reference_data/human_ancestor_GRCh37_e59/human_ancestor_$i.fa.gz | bgzip -c > /net/bipolar/jedidiah/testpipe/vcfs/ancestral/chr$i.aa.vcf.gz";
	&forkExecWait($processcmd);
	
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
    }
    elsif ( $kidpid == 0 ) {
		exec($cmd);
		die "Cannot exec $cmd: $!";
    }
    else {
		waitpid($kidpid,0);
    }
}