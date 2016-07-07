#!/bin/perl

##############################################################################
# Macro for quickly adding, committing, and pushing changes to git repo
# Must be run from within root folder of repo
#
# To-do:
# [-] enable commit message and branch to be modified as options
# [-]
##############################################################################

use strict;
use warnings;
use POSIX;
use Getopt::Long;
use Pod::Usage;

my $branch="master";
my $m="process glfs";
my $help=0;
my $man=0;

# Set options and inputs
GetOptions ('m=s'=> \$m,
'b=s' => \$branch,
'help|?'=> \$help,
man => \$man) or pod2usage(1);

my $datestring = localtime();
my $commit="$datestring: $m";

my $addcmd=`git add *`;
my $commitcmd=`git commit -m \"$commit\"`;
my $pushcmd=`git push origin $branch`;
