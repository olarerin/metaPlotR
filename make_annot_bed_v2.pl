#!/usr/bin/perl
use strict;
use warnings;
use Bio::Perl;
use Bio::SeqIO;
use Getopt::Long;

# Store command line arguements
 my $genomeDir = '';
 my $genePred = '';	
 GetOptions ('genomeDir' => \$genomeDir, 'genePred' => \$genePred); #maybe add option for single multifasta file

# Check command line args are minimally valid
&checkArgs($genomeDir, $genePred) or die;

# Annotate every genic nucleotide
open (ANNOT, "$genePred")  || die "Cant open file: $!";

my $counter = 0;
my $header = <ANNOT>;
my $chrom_seq;
my $last_chr = '';

while (my $line = <ANNOT>){
	chomp ($line);

	# Split line and assign variables
	my @temp = split(/\t/, $line);
	my ($transcript_ID, $chr, $strand) = @temp[1..3];
	my ($txStart, $txEnd, $cdsStart, $cdsEnd) = @temp[4..7];
	my $exonCount = $temp[8];
	my @exonFrame = split(/,/, $temp[-1]);
	my @exon_starts = split (/,/, $temp[9]);
	my @exon_ends = split (/,/, $temp[10]);
	my $gene_name = $temp[12];

	# Load genomic sequence of chromosome (if necessary)
	if ($chr ne $last_chr){
		$chrom_seq = loadGenome($genDir, $chr);
		$last_chr = $chr;
	}

	# Determine size of transcript
	my $gene_size = 0;
	my $num_exons = scalar (@exon_starts);
	for (my $i = 0; $i < $num_exons; $i++){
		$gene_size += abs ($exon_ends[$i] - $exon_starts[$i]);
	}

	# Conversion from 0-based to 1-based coordinates
	$txStart = $txStart + 1;
	$cdsStart = $cdsStart + 1;
	foreach my $x (@exon_starts){
		$x = $x + 1
	}

	# Label the exonic loci
	my $nucl_j = -1;
	my $count = 0;
	my @temp_out = ();

	for (my $i = 0; $i < $exonCount; $i++){ # Go through each exon delimiter
		foreach my $coord (($exon_starts[$i])..($exon_ends[$i])){ # Go through each position of the exon
			$count++;
			my $mrna_pos = &determine_mrna_pos($strand, \$gene_size, $count);
			my $feature = &utrs_or_cds($coord, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd);
			my $nucleotide = uc(substr ($chrom_seq, $coord-1, 1));
			my $codon_pos = &determine_codon_pos($feature, $nucl_j, $strand);

			print $chr, "\t", $coord-1, "\t", $coord, "\t", join('|', ($transcript_ID, $gene_name, $feaure, $mrna_pos)), "\t";
			print $nucleotide, "\t", $strand, "\n";
		}
	}
}

sub determine_codon_pos(){
	my ($feature, $nucl_j, $strand) = @_;
	my $codon_pos = 'NA';

	if ($feature eq 'cds'){
		$nucl_j++;
		$codon_pos = $nucl_j % 3;

		# correct $codon_pos for strandedness
		if ($strand eq '-'){
			if ($codon_pos == 0){ $codon_pos = 2; }
			elsif ($codon_pos == 2){ $codon_pos = 0; }
		}
	}
	return $codon_pos;
}

sub determine_mrna_pos(){
	my ($strand, $gene_size_ref, $count) = @_;
	my $mrna_pos = $count;
	if ($strand eq '-'){
		$mrna_pos = $$gene_size_ref;
		$$gene_size_ref -= 1;
	}  
	return $mrna_pos;
}

sub checkArgs(){
	my ($genDir, $genPred) = @_;
	my pass;
	if (! -d $genDir){
		print "No such directory exists: $genDir \n";
		$pass = 0;
	}
	if (! -f $genPred){
		print "No such file exists: $genDir \n";
		$pass = 0
	}
	return ($pass);
}

sub loadGenome(){
	my ($genDir, $chr) = @_;
	my $file = $genDir.'/'.$chr.'.fa';
	my $seqio_obj = Bio::SeqIO->new(-file => $dir.$file);
	my $seq_obj = $seqio_obj->next_seq;
	return ($seq_obj->seq);
}

sub utrs_or_cds (){
	my ($coord, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd) = @_;
	my $feature = 'na';
	my ($utr5, $cds, $utr3) = (0,0,0);

	if ($coord >= $txStart && $coord < $cdsStart){
		if ($strand eq '+'){$feature = '5utr';}
		else {$feature = '3utr';}
	}
	elsif ($coord > $cdsEnd && $coord <= $txEnd){
		if ($strand eq '+'){$feature = '3utr';}
		else {$feature = '5utr';}
	}
	else{$feature = 'cds';}

	return ($feature);
}
