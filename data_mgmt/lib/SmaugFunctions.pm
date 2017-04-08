package SmaugFunctions;

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(forkExecWait getRef);

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
sub getRef{
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
# read array of filenames and returns file handles
##############################################################################
sub get_write_handles {
  my @file_names = @_;
  my %file_handles;
  foreach (@file_names) {
    open my $fh, '>', $_ or next;
    $file_handles{$_} = $fh;
  }
  return %file_handles;
}

1;
