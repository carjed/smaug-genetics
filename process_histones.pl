#!/usr/local/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Cwd;


# my $bcftools="/net/bipolar/jedidiah/bcftools/bcftools";
# my $vcftools="/net/bipolar/jedidiah/vcftools_0.1.10/bin/vcftools";


my @bedfiles=</net/bipolar/jedidiah/mutation/reference_data/histone_marks/*.bed>;
my @bedindex;
foreach my $file (@bedfiles){
	my $filename=basename($file, ".bed");
	my @strs=split(/\./, $filename);
	my $name= $strs[2];
	print "Binning $name ...\n";
	# print "$filename\n";
	
	my $cmd= "/net/bipolar/jedidiah/bin/bedmap --echo --count --fraction-map 0.51 /net/bipolar/jedidiah/mutation/reference_data/genome.100kb.bed $file > /net/bipolar/jedidiah/mutation/reference_data/histone_marks/binned/$name.100kb.bed";
	
	&forkExecWait($cmd);
	
	print "Done\n";
	
}


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