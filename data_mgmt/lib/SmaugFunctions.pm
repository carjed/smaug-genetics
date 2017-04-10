package SmaugFunctions;

use strict;
use warnings;
use Exporter qw(import);
no warnings 'experimental::smartmatch';
our @EXPORT_OK = qw(forkExecWait getRef getMotif getType getWriteHandles);

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
# Get reference
##############################################################################
sub getRef {
	my $f_fasta=shift;
  my $chr=shift;
  my $nextchr;

  if ($chr<22) {
  	$nextchr=$chr+1;
  } elsif ($chr==22) {
  	$nextchr="X";
  } else {
  	$nextchr="Y";
  }

	open my $fasta, '<', $f_fasta or die "can't open $f_fasta: $!";

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

##############################################################################
# Get K-mer motif
##############################################################################
sub getMotif {
  my $localseq=shift;
  my $adj=shift;
  my $subseq=$adj*2+1;

  my $altlocalseq = reverse $localseq;
  $altlocalseq  =~ tr/ACGT/TGCA/;

  my $ref1 = substr($localseq, $adj, 1);
  my $ref2 = substr($altlocalseq, $adj, 1);

  my $seqp;
  if($ref1 ~~ [qw( A C )]){
    $seqp = "$localseq\($altlocalseq\)";
  } else {
    $seqp = "$altlocalseq\($localseq\)";
  }

  return $seqp;
}

sub getType {
  my $ref=shift;
  my $alt=shift;
  my $adj=shift;
  my $seqp=shift;

  my $CAT = "${ref}${alt}";
  my $Category;
  if($CAT ~~ [qw( AC TG )]){ $Category = "AT_CG";}
  elsif($CAT ~~ [qw( AG TC )]){ $Category = "AT_GC";}
  elsif($CAT ~~ [qw( AT TA )]){ $Category = "AT_TA";}
  elsif($CAT ~~ [qw( GA CT )]){ $Category = "GC_AT";}
  elsif($CAT ~~ [qw( GC CG )]){ $Category = "GC_CG";}
  elsif($CAT ~~ [qw( GT CA )]){ $Category = "GC_TA";}

  if(substr($seqp, $adj, 2) eq "CG"){ $Category = "cpg_$Category";}
  return $Category;
}

##############################################################################
# read array of filenames and returns file handles
##############################################################################
sub getWriteHandles {
  my @file_names = @_;
  my %file_handles;
  foreach (@file_names) {
    open my $fh, '>', $_ or next;
    $file_handles{$_} = $fh;
  }
  return %file_handles;
}

1;
