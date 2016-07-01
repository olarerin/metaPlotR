#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Store command line arguements
 my $annot = '';	
 GetOptions ('annot=s' => \$annot); #maybe add option for single multifasta file

# Check command line args are minimally valid
&checkArgs($annot) || die "Check command line arguments: $!";

# Determines the size of all 5'UTR, CDS, and 3'UTRs in the annotation file
open (ANNOT, "$annot") || die "Cant open file: $!";
#open (OUT, ">/home/shareowner/annotations/d_melanogaster/dme3_UTR_cds_starts_stops.txt") || die "Cant open file: $!";

my %region_start_stop_h;
my %geneID_gname_h;
my %trx_len_h;

# This loop iterates over the annotation bedfile generated from "make_annot_bed.pl" and determines the stop and start sites
# for the 5UTR, CDS and 3UTR in transcriptomic space and saves this in the $region_start_stop hash

while (<ANNOT>){
	chomp;
	my $trx_desc = ( split(/\t/, $_) )[3];
	my ($geneID, $gname, $region, $mRNA_pos) = split (/\|/, $trx_desc);
	$geneID_gname_h{ $geneID } = $gname;

	if ( !exists $trx_len_h{$geneID} ){
		$trx_len_h { $geneID } = 0;
	}
	$trx_len_h{ $geneID }++;

	if ( exists $region_start_stop_h{$geneID}{$region} ){ 	
		my @start_stop = @{$region_start_stop_h{$geneID}{$region}}; # $start_stop[0] = minimum, $start_stop[1] = maximum
		if ($start_stop[0] > $mRNA_pos){
			$start_stop[0] = $mRNA_pos; 
		}
		if ($start_stop[1] < $mRNA_pos){
			$start_stop[1] = $mRNA_pos; 
		}
		$region_start_stop_h{$geneID}{$region} = [@start_stop];
	}
	else{
		$region_start_stop_h{$geneID}{$region} = [($mRNA_pos, $mRNA_pos)];
	}
}


foreach my $geneID ( keys %region_start_stop_h ){
	print $geneID, "\t", $geneID_gname_h{ $geneID }, "\t", $trx_len_h{ $geneID }, "\t";
	if (exists $region_start_stop_h{ $geneID }{ '5utr' }){
		print join ("\t", @{ $region_start_stop_h{ $geneID }{ '5utr' }}), "\t";
	}
	else{
		print "NA\tNA\t";
	}

	if (exists $region_start_stop_h{ $geneID }{ 'cds' }){
		print join ("\t", @{ $region_start_stop_h{ $geneID }{ 'cds' } }), "\t";
	}
	else{
		print "NA\tNA\t";
	} 

	if (exists $region_start_stop_h{ $geneID }{ '3utr' }){
		print join ("\t", @{ $region_start_stop_h{ $geneID }{ '3utr' } }), "\n";
	}
	else{
		print "NA\tNA\n";
	}
}

######################### sub routunes ##########################
sub checkArgs(){
	my ($genPred) = @_;
	my $pass = 1;
	if (! -f $genPred){
		print "No such file exists: $genPred \n";
		$pass = 0;
	}
	return ($pass);
}
