#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Determines the relative and absolute distances of the loci in the annotated
# bed file from the start and end coordinates of the 5'UTR, cds and 3'UTR

# Store command line arguements
 my $bed;
 my $regions;	
 my $abs_out;
 my $rel_out;

GetOptions ('bed=s' => \$bed, 'regions=s' => \$regions);

open (UTRCDS, "$regions") || die "Can't open cds/utr definition file: $!";
open (BED, "$bed") || die "Can't open bed file: $!";

my %refseq_endpts_h;

# Save array of endpoints in hash
# print 'Storing coordinates for UTR and CDS starts and ends..', "\n";
while (my $line = <UTRCDS>){
	next if $line =~ m/(\tNA){1}/;  # IMPORTANT: excludes ncRNA (they lack 5'UTR, cds, and 3'UTR annotations)
	chomp($line);
	my @temp = split (/\t/, $line);
	my $refseq = $temp[0];
	my $gene_symbol = $temp[1];
	my $trx_len = $temp[2];
	$refseq_endpts_h{$refseq} = [@temp[3..$#temp]];
}

# Iterate through annotated bed file and determine the relative
# and absolute distance to the endpoints

# Print headers for output files
print "chr\tcoord\tgene_name\trefseqID\trel_location\t";
print "utr5_st\tutr5_end\tcds_st\tcds_end\tutr3_st\tutr3_end\t";
print "utr5_size\tcds_size\tutr3_size\n";

#print 'Iterating through bed file...', "\n";

my $excluded = 0;
my $total = 0;

while (my $line = <BED>){
	next if $line =~ m/^\#/;
	chomp ($line);
	$total++;
	#my ($chr, $end, $gname, $refseq, $mrna_pos) = (split(/\t/, $line))[0,2,3,6,7];
	my ($chr, $end, $mrna_meta) = (split(/\t/, $line))[6,8,9];
	my ($refseq, $gname, $region, $mrna_pos) = split (/\|/, $mrna_meta);

	if (!exists $refseq_endpts_h{$refseq}){
		$excluded++;
		next;
	}
	my @endpts = @{$refseq_endpts_h{$refseq}};

	my @abs_dist = &abs_distance($mrna_pos, \@endpts); 
	my $rel_dist = &rel_distance($mrna_pos, \@endpts); 
	my @feature_sizes = &size_features(\@endpts);

	print "$chr\t$end\t$gname\t$refseq\t", $rel_dist, "\t", join("\t", @abs_dist), "\t";
	print join("\t", @feature_sizes), "\n";

}

print STDERR '** Total bedfile lines: ', "$total\n";
print STDERR '** Excluded (features not well defined or gene not found in gtf): ', "$excluded\n";

######################### sub routunes ##########################
sub abs_distance {
	my ($mrna_pos, $endpts_ref) = @_;
	my @output;

	foreach my $point (@{$endpts_ref}){
		push @output, ($mrna_pos - $point);
	}
	return @output;
}

sub rel_distance {
	my ($mrna_pos, $endpts_ref) = @_;
	my ($u5_st, $u5_end, $cds_st, $cds_end, $u3_st, $u3_end) = @{$endpts_ref};
	my @output;
	my $rel_pos;

	if ($mrna_pos >= $u5_st and $mrna_pos <= $u5_end){
		$rel_pos = ( $mrna_pos-$u5_st + 1) / ( $u5_end-$u5_st + 1); 
	}
	elsif ($mrna_pos >= $cds_st and $mrna_pos <= $cds_end){
		$rel_pos = ($mrna_pos - $cds_st + 1)/($cds_end - $cds_st + 1); 
		$rel_pos += 1;
	}
	elsif ($mrna_pos >= $u3_st and $mrna_pos <= $u3_end){
		$rel_pos = ( $mrna_pos-$u3_st + 1) / ( $u3_end-$u3_st + 1); 
		$rel_pos += 2;
	}
	return $rel_pos;
}

sub size_features{
	my ($endpts_ref) = @_;
	return ($$endpts_ref[1]-$$endpts_ref[0],
		$$endpts_ref[3]-$$endpts_ref[2],
		$$endpts_ref[5]-$$endpts_ref[4]); 
}


