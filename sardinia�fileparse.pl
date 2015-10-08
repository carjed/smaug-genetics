#!/usr/local/bin/perl

##############################################################################
# Development script for parsing filenames in perl
##############################################################################

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Cwd;



my @vcfs = </net/bipolar/lockeae/freeze4/vcfs/anno/final/*.vcf.gz>;



	foreach my $file (@vcfs) {
		my $filename=fileparse($file);
		my $subfile = substr($filename, index($filename, 'chr'), index($filename, 'anno'));
		my $chr = substr($subfile, 0, index($subfile, '.'));
		print "$chr\n";
		#my $tabix ="tabix -r header.txt in.vcf.gz > out.vcf.gz";
		#&forkExecWait($tabix);
	}
