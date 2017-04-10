#!/usr/bin/env perl

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
use SmaugFunctions qw(getMotif getType);

my $opts = parse_params();
fill_info($opts,$$opts{aa_file});

exit;

#--------------------------------

sub error {
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    die
      "About: This script fills sequence motif and mutation type into the INFO column of VCF files \n";
}


sub parse_params {
    my $opts = {};
    while (my $arg=shift(@ARGV)) {
        if ( $arg eq '-a' || $arg eq '--ancestral-allele' ) { $$opts{aa_file} = shift(@ARGV); next }
        if ( $arg eq '-t' || $arg eq '--type' ) {
            my %known = ( snp=>'s', indel=>'i', all=>'a', ref=>'r' );
            my $types = shift(@ARGV);
            for my $t (split(/,/,$types)) {
                if ( !(exists($known{$t})) ) { error("Unknown type [$t] with -t [$types]\n"); }
                $$opts{types}{$known{$t}} = 1;
            }
            if ( exists($$opts{types}{a}) ) {
                $$opts{types}{s} = 1;
                $$opts{types}{i} = 1;
                $$opts{types}{r} = 1;
            }
            next;
        }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( !exists($$opts{aa_file}) ) { error("Missing the -a option.\n") }
    return $opts;
}


sub fill_info {
    my ($opts,$fa_fname) = @_;

    my $n_unknown = 0;
    my $n_filled_sites = 0;
    my $n_filled_bases = 0;

    my $vcf = Vcf->new(fh=>\*STDIN, assume_uppercase=>1);
    $vcf->parse_header();
    $vcf->add_header_line({key=>'INFO',ID=>'Motif',Number=>1,Type=>'String',
        Description=>'K-mer sequence motif centered at site'});
    print $vcf->format_header();

    my %chr2fa = ();
    my $nskipped = 0;
    while (my $line = $vcf->next_line() )
    {
        my $rec = $vcf->next_data_array($line);
        my $chr = $$rec[0];
        my $pos = $$rec[1];
        my $ref = $$rec[3];
        my $alt = $$rec[4];

        if (!exists($chr2fa{$chr})) {
            my $fname = $fa_fname;
            if ( ! -e $fname ) {
                if ( -e "$fname$chr.fasta.gz" ) { $fname = "$fname$chr.fasta.gz"; }
                else { error(qq[Neither "$fname" nor "$fname$chr.fasta.gz" exists.\n]); }
            }
            $chr2fa{$chr} = FaSlice->new(file=>$fname, size=>100_000);
        }

        my $fa = $chr2fa{$chr};
        my $ref_len = length($ref);
        if ( exists($$opts{types}) && !exists($$opts{types}{a}) ) {
            my $ok = 0;
            for my $alt (split(/,/,$$rec[4])) {
                my ($type,$len,$ht) = $vcf->event_type($ref,$alt);
                if ( exists($$opts{types}{$type}) ) { $ok=1; last; }
            }
            if ( !$ok ) {
                print $line;
                $nskipped++;
                next;
            }
        }

        my $motif = $fa->get_slice($chr, $pos-$adj, $pos+$adj);
        $motif = getMotif($motif, $adj);
        if ( $motif ) {
            $$rec[7] = $vcf->add_info_field($$rec[7],'Motif'=>$motif);
            $n_filled_sites++;
            $n_filled_bases+=$ref_len;
        } else {
            $$rec[7] = $vcf->add_info_field($$rec[7],'Motif'=>'.');
            $n_unknown++;
        }

        my $type = getType($ref, $alt, $adj, $motif);
        # $aa = getMotif($aa, $adj);
        if ( $type ){
            $$rec[7] = $vcf->add_info_field($$rec[7],'Category'=>$type);
        } else {
            $$rec[7] = $vcf->add_info_field($$rec[7],'Category'=>'.');
        }

        print join("\t",@$rec),"\n";
    }

    print STDERR
        "Motif sites filled  .. $n_filled_sites\n",
        "Motif bases filled  .. $n_filled_bases\n",
        "No motifs           .. $n_unknown\n",
        "Lines skipped    .. $nskipped\n";
}
