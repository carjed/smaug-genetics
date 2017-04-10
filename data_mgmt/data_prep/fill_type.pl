#!/usr/bin/env perl
#
# Notes:
#   * The AA files can be downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/technical/reference/ancestral_alignments
#   * The program runs samtools, therefore the AA files must be gzipped (not b2zipped).
#
# support: pd3@sanger

use strict;
use warnings;
use Carp;
use Vcf;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin";
use FaSlice;
use YAML::XS 'LoadFile';
use feature 'say';

my $relpath = $FindBin::Bin;
my $configpath = dirname(dirname($relpath));
my $config = LoadFile("$configpath/_config.yaml");

# print "Script will run with the following parameters:\n";
# for (sort keys %{$config}) {
#     say "$_: $config->{$_}";
# }

my $adj = $config->{adj};
my $mac = $config->{mac};
my $binw = $config->{binw};
my $data = $config->{data};
my $bin_scheme = $config->{bin_scheme};
my $parentdir = $config->{parentdir};
my $count_motifs = $config->{count_motifs};
my $expand_summ = $config->{expand_summ};

use lib "$FindBin::Bin/../lib";
use SmaugFunctions qw(forkExecWait getRef getMotif getType);

# my $opts = parse_params();
# fill_type($opts,$$opts{aa_file});
fill_type();

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    die
        "About: This script fills ancestral alleles into INFO column of VCF files. It depends on samtools,\n",
        "   therefore the fasta sequence must be gzipped (not bgzipped!) and indexed by samtools faidx.\n",
        "   The AA files can be downloaded from\n",
        "       ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/technical/reference/ancestral_alignments\n",
        "   and processed as shown in the example below. This is because the sequences in the original files\n",
        "   are named as 'ANCESTOR_for_chromosome:NCBI36:1:1:247249719', but the underlying FaSplice.pm\n",
        "   requires names as 'chr1' or '1'.\n",
        "Usage: fill-aa [OPTIONS] < in.vcf >out.vcf\n",
        "Options:\n",
        "   -a, --ancestral-allele <prefix>     Prefix to ancestral allele chromosome files.\n",
        "   -t, --type <list>                   Variant types to process: all,indel,ref,snp. [all]\n",
        "   -h, -?, --help                      This help message.\n",
        "Example:\n",
        "   # Get the files ready: compress by gzip and index by samtools faidx. Either repeat the\n",
        "   # following command for each file manually\n",
        "   bzcat human_ancestor_1.fa.bz2 | sed 's,^>.*,>1,' | gzip -c > human_ancestor_1.fa.gz\n",
        "   samtools faidx human_ancestor_1.fa.gz\n",
        "   \n",
        "   # .. or use this loop (tested in bash shell)\n",
        "   ls human_ancestor_*.fa.bz2 | while read IN; do\n",
        "       OUT=`echo \$IN | sed 's,bz2\$,gz,'`\n",
        "       CHR=`echo \$IN | sed 's,human_ancestor_,, ; s,.fa.bz2,,'`\n",
        "       bzcat \$IN | sed \"s,^>.*,>\$CHR,\" | gzip -c > \$OUT\n",
        "       samtools faidx \$OUT\n",
        "   done\n",
        "   \n",
        "   # After this has been done, the following command should return 'TACGTGGcTGCTCTCACACAT'\n",
        "   samtools faidx human_ancestor_1.fa.gz 1:1000000-1000020\n",
        "   \n",
        "   # Now the files are ready to use with fill-aa. Note that the VCF file\n",
        "   # should be sorted (see vcf-sort), otherwise the performance would be seriously\n",
        "   # affected.\n",
        "   cat file.vcf | fill-aa -a human_ancestor_ 2>test.err | gzip -c >out.vcf.gz \n",
        "\n";
}

sub fill_type
{
    # my ($opts,$aa_fname) = @_;

    my $n_unknown = 0;
    my $n_filled_sites = 0;

    my $vcf = Vcf->new(fh=>\*STDIN, assume_uppercase=>1);
    $vcf->parse_header();
    $vcf->add_header_line({key=>'INFO',ID=>'Category2',Number=>1,Type=>'String',
        Description=>'Mutation type'});
    print $vcf->format_header();

    my %chr2fa = ();
    while (my $line = $vcf->next_line() )
    {
        my $rec = $vcf->next_data_array($line);
        my $chr = $$rec[0];
        my $pos = $$rec[1];
        my $ref = $$rec[3];
        my $alt = $$rec[4];

        my $type = getType($ref, $alt);
        # $aa = getMotif($aa, $adj);
        if ( $type )
        {
            $$rec[7] = $vcf->add_info_field($$rec[7],'Category2'=>$type);
            $n_filled_sites++;
        }
        else
        {
            $$rec[7] = $vcf->add_info_field($$rec[7],'Category2'=>'.');
            $n_unknown++;
        }
        print join("\t",@$rec),"\n";
    }

    print STDERR
        "Types filled  .. $n_filled_sites\n",
        "No data           .. $n_unknown\n",
}
