#!/bin/perl

##############################################################################
# Master script for glf processing
##############################################################################

use strict;
use warnings;
use POSIX;

my $branch="develop";

my $datestring = localtime();
# print "Batch job $ID queued at $datestring...\n";
my $commit="$datestring: process glfs";

my $addcmd=`git add *`;
my $commitcmd=`git commit -m \"$commit\"`;
my $pushcmd=`git push origin $branch`;
