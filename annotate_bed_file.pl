#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#  This script is just a wrapper to run intersectBed on your bedfile to annotate each nucleotide position

# Store command line arguements
 my $bed;
 my $annot;	
 GetOptions ('bed=s' => \$bed, 'bed2=s' => \$annot);

my $cmd = 'sort -k1,1 -k2,2n ';
$cmd .= $bed;
$cmd .= '  | intersectBed -a stdin -b ';
$cmd .= $annot;
$cmd .= ' -sorted -wo -s ';

print STDERR 'Run command: ', $cmd, "\n";
system($cmd);
