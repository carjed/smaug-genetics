#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Temp qw( :POSIX );

use vars qw($VERSION);
$VERSION = '0.03';

# Try to load Pod::Usage and install a fallback if it doesn't exist
eval {
    require Pod::Usage;
    Pod::Usage->import();
    1;
} or do {
    *pod2usage = sub {
        die "Error in command line.\n";
    };
};

GetOptions(
    "disk" => \my $do_tie,
    "on|j=s" => \my $joincol,
    "left|j1|1=s" => \my @left_key_cols,
    "right|j2|2=s" => \my @right_key_cols,
    "output|o" => \my @output_fieldlist,
    "delimiter|t|d=s" => \my $delimiter,
    "output-delimiter|od" => \my $output_delimiter,
    "missing|v=i" => \my @missing,
    "null|n=s" => \my $nullvalue,
    "warn-on-duplicates|u=s" => \my @warn_on_duplicates,
    "die-on-duplicates|s=s" => \my @die_on_duplicates,
    "smart-duplicates" => \my $smart_duplicates,
    "progress|verbose|p" => \my $progress,
    'help'              => \my $help,
    'version'           => \my $version,
) or pod2usage(2);
pod2usage(1) if $help;
if (defined $version) {
    print "$VERSION\n";
    exit 0;
};
pod2usage("$0: No files given.")  if ((@ARGV == 0) && (-t STDIN));

$delimiter ||= "\t";
$nullvalue ||= "";
$output_delimiter ||= $delimiter;
$joincol ||= 1;
if (! @left_key_cols) { @left_key_cols = $joincol; };
if (! @right_key_cols) { @right_key_cols = $joincol; };
my %output_missing = map { $_ => $_ } @missing;

my %col_count;

for (\@left_key_cols, \@right_key_cols) {
    @$_ = map { split /,/ } @$_;
};
if (@left_key_cols != @right_key_cols) {
    local $" = ",";
    warn "Left  keys: @left_key_cols\n";
    warn "Right keys: @right_key_cols\n";
    die "Differing number of key columns between left and right - that is wrong.\n";
};

# Adjust the indices for the join columns:
for (@left_key_cols, @right_key_cols) {
    $_--
};

my @output_cols = ('1.*','2.%');
if (@output_fieldlist) {
    @output_cols = map { split /,/ } @output_fieldlist;
};

if (@warn_on_duplicates and not @die_on_duplicates) {
    @warn_on_duplicates = (2);
};

my %on_duplicate;
$on_duplicate{ $_ } = sub { warn "Duplicate key '$_[0]' for row >>$_[1]<< in file $_[2]\n" } for @warn_on_duplicates;
$on_duplicate{ $_ } = sub { die  "Duplicate key '$_[0]' for row >>$_[1]<< in file $_[2]\n" } for @die_on_duplicates;

my %right; # The index into the right file
my %seen;  # The keys we processed from the left file

my @CLEANUP;
if ($do_tie) {
    require DB_File;
    my $rn = tmpnam;
    tie %right, 'DB_File', $rn;
    my $sn = tmpnam;
    tie %seen, 'DB_File', $sn;
    push @CLEANUP, $rn, $sn;
};
END {
    if ($do_tie) {
        untie %right;
        untie %seen;
    };
    for (@CLEANUP) {
        unlink $_ or warn "Couldn't remove tempfile '$_' : $!\n";
    };
};

my ($left,$right) = @ARGV;

# Read the right file into the hash
open my $rfh, "<", $right
    or die "Couldn't read '$right': $!";
open my $lfh, "<", $left
    or die "Couldn't read '$left': $!";

sub key {
    my ($cols,$col_info) = @_;
    return join $delimiter, @{ $cols }[ @$col_info ];
};

sub output {
    my @lr = (@_);
    my @output = map { /^(\d)\.(\d+)/ or die "Invalid column spec '$_'"; $lr[$1-1]->[$2-1] } @output_cols;
    print join($delimiter, @output), "\n";
};

sub expand_output_columns {
    my (@list) = @_;

    my %keycols = (
         1 => +{ map { $_+1 => 1 } @left_key_cols },
         2 => +{ map { $_+1 => 1 } @right_key_cols },
    );

    my @res = map { /(\d)\.\*/ ? (map { "$1.$_" } (1..$col_count{ $1 }))
                  : /(\d)\.\%/ ? (map { "$1.$_" } grep { ! exists $keycols{$1}{$_}} (1..$col_count{ $1 }))
                  : $_
              } @list;
    @res
};

warn "Reading $right"
    if $progress;
while (<$rfh>) {
    chomp;
    my @right_cols = split /\Q$delimiter\E/;
    $col_count{ 2 } ||= @right_cols;
    my $key = key( \@right_cols, \@right_key_cols );
    if ($right{ $key } and $on_duplicate{2}) {
        my $diff = $right{ $key } ne $_;
        if ($diff or !$smart_duplicates) {
            $on_duplicate{2}->($key,$_,$right)
        };
    };
    $right{ $key } = $_;
};

# Read the left file and output the generated lines (if any)
warn "Processing $left"
    if $progress;
my $expanded_output_columns;
while (<$lfh>) {
    chomp;
    my @left_cols = split /\Q$delimiter\E/;
    $col_count{ 1 } ||= @left_cols;
    my $key = key( \@left_cols, \@left_key_cols );

    if ($seen{ $key } and $on_duplicate{1}) {
        my $diff = $seen{ $key } ne $_;
        if ($diff or !$smart_duplicates) {
            $on_duplicate{1}->($key,$_,$left)
        };
    };
    $seen{ $key }++;
    my $out;
    my @right_cols;
    if (exists $right{ $key }) {
        @right_cols = split /\Q$delimiter\E/, $right{ $key };
    } else {
        @right_cols = ($nullvalue) x $col_count{ 2 };
    };

    if (exists $right{ $key } or $output_missing{1}) {
        if (! $expanded_output_columns) {
            @output_cols = expand_output_columns(@output_cols);
            $expanded_output_columns++;
        };
        output \@left_cols, \@right_cols;
    };
};

@output_cols = expand_output_columns(@output_cols);

if ($output_missing{2}) {
    warn "Writing right-missing keys"
        if $progress;
    my @left_cols = ($nullvalue) x $col_count{ 1 };
    while ((my ($key,$v)) = each %right) {
        if (! $seen{ $key }) {
            my @right_cols = split /\Q$delimiter\E/, $v;
            output \@left_cols, \@right_cols;
        };
    };
};

__END__

=head1 NAME

join - join two files by common key columns

=head1 SYNOPSIS

  join.pl [OPTIONS] FILE1 FILE2

  join.pl --on 1,2 file1.txt file2.txt

  join.pl --left 1,2 --right 3,4 file1.txt file2.txt

=head1 OPTIONS

=item B<--on COL> - specify a single column number to join both files on

This is a shorthand for C<--left COL --right COL>

=item B<--missing FILE> - output rows only in one file

C<--missing 1> will output rows that only exist in the left file.

=item B<--null VAL> - string for the null value

When a row is output through the C<--missing> option, the missing
values will be replaced by the value given.

The default is an empty string, "".

Example: --null NULL

=item B<--warn-on-duplicates FILE> - output a warning if duplicate keys are found in the file

=item B<--die-on-duplicates FILE> - die if duplicate keys are found in the file

These options govern how the program behaves when it encounters
duplicate keys in a file.

=item B<--smart-duplicates> - be smart about duplicates

This setting enables smart duplicate handling that will
only consider a row as duplicate if the key is identical but
the remaining values differ.

=item B<--left COL1,COL2> - specify key columns for the left file

=item B<--right COL1,COL2> - specify key columns for the right file

The column counts starting at 1. The default column is 1.

=item B<--output COL1,COL2> - specify columns to output

If you want to reorder or omit columns use this to
list the columns. Each column must be in the format
C<M.N> where M is either 1 for the left file or
2 for the right file, and N is the column number.

There are two shorthands:

C<M.*> will include all columns from the source file
in source order.

C<M.%> will include all columns from the source file
except the key columns in source order.

The default is C<1.* 2.%>, which will append
the non-key columns of the right file to the left file.

=item B<--delimiter DEL> - specify column input delimiter

The default input delimiter is a tab. No automatic
delimiter recognition is done yet.

=item B<--output-delimiter DEL> - specify output column delimiter

The output column delimiter defaults to the input
column delimiter.

=item B<--progress> - be verbose in the progress

Some diagnostic messages will be output to STDERR
as the program progresses.

=item B<--disk> - use disk memory for joining instead of RAM

This will use disk memory for storing the index instead
of using RAM.

=item B<--version> - print program version

Outputs the program version.

=item B<--help> - print this page

Outputs this help text.

