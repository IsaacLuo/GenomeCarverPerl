#! /usr/bin/perl
#2012-02-15
#v2 changes:
#1. Only verified genes are calculated
#2. Intron sequences are marked down with lower case
#3. Introns are calculated from GFF file, but not from the UCSC file
#4. Essential genes are from SGD YeastMine, not from BioStudio

#2012-03-08
#bug fixed: 
#Reverse gene intron sequnece is now corrected
#3UTR length is 200 but not 199
#Promoter length is 500 but not 499


#2014-1-16
#Re-written the code for designing oligos for yeastfab project

use strict;
use warnings;

###Change the file name here!!
my $file_name = "YeastFab";


#my $query_file = shift; 


my $gff_file = "/Users/Cai/Dropbox/Programming/BitBucket/Genome_Carver/Yeast_chromosomes/saccharomyces_cerevisiae.gff";

#open(FILE, "<", $query_file) or die "Cannot open file $query_file!\n";

#my @query_genes = <FILE>;

#open(INTRON, "<", $intron_file) or die "Cannot open file $intron_file!\n";
#
#my @introns = <INTRON>; 

#my @query_genes=();
#
#while (my $gene = <FILE>) {
#	push(@query_genes, $gene);
#}

open(GFF, "<", $gff_file) or die "Cannot open file $gff_file!\n";

#Get the Chromosome sequences

my $entire_gff;

while (my $line = <GFF>) {
	$entire_gff .= $line;
}

my @entire_gff = split(/\#\#FASTA/, $entire_gff); 

my $fasta_section = $entire_gff[1];

my @fasta_section = split(/>/, $fasta_section); 

my %chromosome_sequence = ();

foreach my $giant_fasta_line(@fasta_section)
{
	my @fasta_lines = split(/\n/, $giant_fasta_line);
	my $chromosome_id; 
	my $chromosome_seq;
	foreach my $line(@fasta_lines) {
		#chomp $line;
		#if ($line =~ m/chr/) {
		#	$chromosome_id = $line;
		#}
		#else 
		#{
		#$chromosome_seq .= $line;
		#}
		if ($line =~ m/chr/) {
			 $chromosome_id = $line;
			 #print $chromosome_id."\n";
		}
		elsif ($line =~ m/^[A|T|C|G|a|t|c|g]+$/)
		{
			$chromosome_seq .= $line;
		}
	}
	#print $chromosome_seq."\n";
	if (defined($chromosome_id))
	{
		$chromosome_sequence{$chromosome_id}->{Seq} = $chromosome_seq;
	}
}

close(GFF);

#up to here, all the chromosomes have their own sequences associated. Good. 

#print $chromosome_sequence{chrI}->{Seq}."\n";



#To get a gene sequence from position X to position Y, the format is {substr $chromosome_sequence, X-1, Y-X+1}
#This gets more complicated with the strand issue

open(GFF, "<", $gff_file) or die "Cannot open file $gff_file!\n";

my $gene_count = 0; 

my %gff_feature = (); 

my %intron_feature = ();
while (my $line = <GFF>) {
	last if ($line =~ m/^\#\#FASTA/); 
	#
	$line =~ s/\%20/ /;
	if ($line =~ m/^\#/) {
		next;
	}
	#my @line_content = split(/\s+/, $line);
	my @line_content = split(/\t/, $line);
	if ($line_content[2] eq "gene") {
		my $id; 
		my @note_content = split (/;/, $line_content[8]);
		foreach my $note_item(@note_content) {
			if ($note_item =~ m/^ID=/)
			{
				#print $'."\n";
				$id = $';
				$gff_feature{$id}->{ID} = $';
			}
			if ($note_item =~ m/^gene=/)
			{
				$gff_feature{$id}->{Gene} = $';
			}
			if ($note_item =~ m/^Alias=/)
			{
				$gff_feature{$id}->{Alias} = $';
			}
			if ($note_item =~ m/^Note=/)
			{
				$gff_feature{$id}->{Note} = $';
			}	
			if ($note_item =~ m/^orf_classification=/)	
			{
				$gff_feature{$id}->{Verification} = $';
				chomp $gff_feature{$id}->{Verification};
			}
		}
#Putting things in places
		$gff_feature{$id}->{Chr} = $line_content[0];
		$gff_feature{$id}->{Source} = $line_content[1];
		$gff_feature{$id}->{Feature} = $line_content[2];
		$gff_feature{$id}->{Start} = $line_content[3];
		$gff_feature{$id}->{End} = $line_content[4];
		$gff_feature{$id}->{Score} = $line_content[5];
		$gff_feature{$id}->{Strand} = $line_content[6];
		$gff_feature{$id}->{Frame} = $line_content[7];
		$gene_count++
	}
	
	if ($line_content[2] eq "intron")
	{
		my $id;
		#chrII	SGD	intron	4117	4215	.	-	.	Parent=YBL111C_mRNA;Name=YBL111C_intron;orf_classification=Uncharacterized		
		my @note_content = split (/;/, $line_content[8]);
		
		foreach my $note_item(@note_content) {
			if ($note_item =~ m/^Name=/)
			{
				
				$id = $';
				if (exists($intron_feature{$id}->{ID}))
				{
					$id .= $line_content[3];
					#append start poisiton to the intron if it is not the only intron in the gene
				}
				$intron_feature{$id}->{ID} = $id."\n";
			}
		}
	
		$intron_feature{$id}->{Chr} = $line_content[0];
		$intron_feature{$id}->{Source} = $line_content[1];
		$intron_feature{$id}->{Feature} = $line_content[2];
		$intron_feature{$id}->{Start} = $line_content[3];
		$intron_feature{$id}->{End} = $line_content[4];
		$intron_feature{$id}->{Score} = $line_content[5];
		$intron_feature{$id}->{Strand} = $line_content[6];
		$intron_feature{$id}->{Frame} = $line_content[7];
		
		#print $intron_feature{$id}->{Chr}."\n";
		#print $intron_feature{$id}->{Source}."\n";
		#print $intron_feature{$id}->{Feature}."\n";
		#print $intron_feature{$id}->{Start}."\n";
		#print $intron_feature{$id}->{End}."\n";
		#print $intron_feature{$id}->{Score}."\n";
		#print $intron_feature{$id}->{Strand}."\n";
		#print $intron_feature{$id}->{Frame}."\n";
		
	}
	
}


#foreach my $key (keys %gff_feature) {
#	print "$gff_feature{$key}->{ID}\t $gff_feature{$key}->{Start}\t $gff_feature{$key}->{End}\n";
#}


######Iterating the gene list to get sequences
#foreach my $query_gene(@query_genes) {
#	chomp $query_gene;
#	print "$gff_feature{$query_gene}->{ID}\t  $gff_feature{$query_gene}->{Chr}\t $gff_feature{$query_gene}->{Start}\t $gff_feature{$query_gene}->{End}\t $gff_feature{$query_gene}->{Alias}\n";
#	print "$gff_feature{$query_gene}->{Note}\n";
#	
#	if (grep (/^$query_gene$/, @introns)) {
#		print "\nThis gene has intron(s)\n";
#	}
#
#	else {
#		print "\nThis gene doesn't have intron\n";
#	}
#	
#	
#	print "\nThe upstream gene is ".search_last_gene($query_gene)."\n";
#
#	print "\nThe downstream gene is ".search_next_gene($query_gene)."\n";
#	
#	#print my $gene_sequence = substr $chromosome_sequence{$gff_feature{YPL160W}->{Chr}}->{Seq}, $gff_feature{YPL160W}->{Start}-1, $gff_feature{YPL160W}->{End}-$gff_feature{YPL160W}->{Start}+1;
#
#	print "\nThe upstream gene is ".search_last_gene($query_gene)."\n";
#
#	print "\nThe downstream gene is ".search_next_gene($query_gene)."\n";
#
#	print "\nThe promoter sequence has ".check_restriction_sites(get_promoter_sequence($query_gene))."\n".get_promoter_sequence($query_gene)."\n";
#
#	print "\nThe ORF sequence has ".check_restriction_sites(get_orf_sequence($query_gene))."\n".get_orf_sequence("$query_gene")."\n";
#
#	print "\nThe 3'UTR sequence has ".check_restriction_sites(get_3UTR_sequence($query_gene))."\n".get_3UTR_sequence($query_gene)."\n***************************************************\n";
#}
#


##############below is to print a CSV file for better inspection#############
my $CSV_file_name = $file_name.".txt"; ##Change file name here
open (CSV, ">", $CSV_file_name);
print CSV "System name\t Chromosome\t Start\t End\t Alias\t Note\t Intron\t ORF Length\t ORF Sequence\t ORF Restriction Sites\t Promoter Sequence\t Promoter Restriction Sites\t 3'UTR sequence\t 3'UTR Restriction Sites\n"; 

my %forward_primers; 
my %reverse_primers; 



#foreach my $key (keys %gff_feature) {
#	print "$gff_feature{$key}->{ID}\t $gff_feature{$key}->{Start}\t $gff_feature{$key}->{End}\n";
#}


foreach my $key (keys %gff_feature) {
	chomp $key;
	
	design_golden_gate_primers($key); 
	
	print CSV "$gff_feature{$key}->{ID}\t $gff_feature{$key}->{Chr}\t $gff_feature{$key}->{Start}\t $gff_feature{$key}->{End}\t $gff_feature{$key}->{Alias}\t $gff_feature{$key}->{Note}\t"; 
	
	my @introns = check_introns($key); 

		
	#if (grep (/^$query_gene$/, @introns)) {
	#	print CSV "YES\t";
	#}
    #
	#else {
	#	print CSV "NO\t";
	#}
	
	print CSV scalar(@introns)."\t";
	
	my $orf_length = length(get_orf_sequence("$key")); 
	
	my $ORF_sequence = get_orf_sequence("$key")."\n";
	
	if ($gff_feature{$key}->{Strand} eq "+")
	{
		foreach my $intron(@introns)
		{
			$ORF_sequence = substr($ORF_sequence, 0, $intron_feature{$intron}->{Start}-$gff_feature{$key}->{Start}).lc(substr($ORF_sequence, ($intron_feature{$intron}->{Start}-$gff_feature{$key}->{Start}), ($intron_feature{$intron}->{End}-$intron_feature{$intron}->{Start})+1)).substr($ORF_sequence, ($intron_feature{$intron}->{End}+1-$gff_feature{$key}->{Start}), ($gff_feature{$key}->{End}-$intron_feature{$intron}->{End}));
			#print $gff_feature{$query_gene}->{Start}."\n";
			#print $gff_feature{$query_gene}->{End}."\n";
			#print $intron_feature{$intron}->{Start}."\n";
			#print $intron_feature{$intron}->{End}."\n";
		}
	}

	#Crick strand -- intron needs special attention
	else
	{
		foreach my $intron(@introns)
		{
			#$ORF_sequence = substr($ORF_sequence, 0, ($gff_feature{$query_gene}->{End}-$intron_feature{$intron}->{End})).lc(substr($ORF_sequence, ($gff_feature{$query_gene}->{End}-$intron_feature{$intron}->{End}), ($intron_feature{$intron}->{End}-$intron_feature{$intron}->{Start}))).substr($ORF_sequence, $gff_feature{$query_gene}->{End}-$intron_feature{$intron}->{Start}, ($intron_feature{$intron}->{Start}-$gff_feature{$query_gene}->{Start}));
			$ORF_sequence = substr($ORF_sequence, 0, ($gff_feature{$key}->{End}-$intron_feature{$intron}->{End})).lc(substr($ORF_sequence, ($gff_feature{$key}->{End}-$intron_feature{$intron}->{End}), ($intron_feature{$intron}->{End}-$intron_feature{$intron}->{Start})+1)).substr($ORF_sequence, ($gff_feature{$key}->{End}-$intron_feature{$intron}->{Start})+1, ($intron_feature{$intron}->{Start}-$gff_feature{$key}->{Start}));
		
			#print $gff_feature{$query_gene}->{Start}."\n";
			#print $gff_feature{$query_gene}->{End}."\n";
			#print $intron_feature{$intron}->{Start}."\n";
			#print $intron_feature{$intron}->{End}."\n";
		}
	
	}
	
	chomp $ORF_sequence;
	
	print CSV "$orf_length\t".$ORF_sequence."\t".check_restriction_sites(get_orf_sequence("$key"))."\t". get_promoter_sequence("$key")."\t".check_restriction_sites(get_promoter_sequence("$key"))."\t".get_3UTR_sequence("$key")."\t".check_restriction_sites(get_3UTR_sequence("$key"))."\n";
	
}


close (CSV); 


#########################produce forward primers#################
my $forward_primer_file_name = $file_name."_forward_primer.txt"; ##Change file name here
open (FWD, ">", $forward_primer_file_name);


foreach my $key (sort keys %forward_primers) 
{
	print FWD $key."\t".$forward_primers{$key}."\n";
}

close FWD;

#########################produce reverse primers#################
my $reverse_primer_file_name = $file_name."_reverse_primer.txt"; ##Change file name here
open (REV, ">", $reverse_primer_file_name);

foreach my $key (sort keys %reverse_primers) {
	print REV $key."\t".$reverse_primers{$key}."\n";
}

close REV;

#########################Subroutines#############################
sub check_introns {
	#This sub takes the gene name and scane for any introns it contains
	my $id = shift;
	my @introns = (); 
	foreach my $key (sort keys %intron_feature) {
		if (($intron_feature{$key}->{Chr} eq $gff_feature{$id}->{Chr}) && ($intron_feature{$key}->{Start} >= $gff_feature{$id}->{Start} && ($intron_feature{$key}->{End} <= $gff_feature{$id}->{End})) && ($key =~ m/^$id/))
		{
			push (@introns, $key);
		}
	}
	return @introns;
}




sub search_last_gene {
	# body...
	my $id = shift; 
	my $distance = 9999999;
	my $last_gene; 
	foreach my $key (keys %gff_feature){
		if (($gff_feature{$id}->{Chr} eq $gff_feature{$key}->{Chr}) && ($gff_feature{$id}->{Start} > $gff_feature{$key}->{End}) && ($gff_feature{$key}->{Verification} =~ m/^Verified/)){
			#check whether two genes are on the same chromosome, and if the query gene is upstream of the target gene; if so, proceed
			if (($gff_feature{$id}->{Start} - $gff_feature{$key}->{End}) < $distance)
			{
				$distance = $gff_feature{$id}->{Start} - $gff_feature{$key}->{End}; 
				$last_gene = $key;
			}
		}
	}
		return $last_gene; 
	
}

sub search_next_gene {
	# body...
	my $id = shift; 
	my $distance = 9999999;
	my $next_gene; 
	foreach my $key (keys %gff_feature){
		if (($gff_feature{$id}->{Chr} eq $gff_feature{$key}->{Chr}) && ($gff_feature{$id}->{End} < $gff_feature{$key}->{Start}) && ($gff_feature{$key}->{Verification} =~ m/^Verified/)){
			#check whether two genes are on the same chromosome, and if the query gene is downstream of the target gene; if so, proceed
			if (($gff_feature{$key}->{Start} - $gff_feature{$id}->{End}) < $distance)
			{
				$distance = $gff_feature{$key}->{Start} - $gff_feature{$id}->{End}; 
				$next_gene = $key;
			}
		}
	}
	return $next_gene; 
}

sub get_orf_sequence {
	# body...
	# This sub takes a System name (i.e., ID) and returns the DNA sequence of it
	my $id = shift;
	if ($gff_feature{$id}->{Strand} eq "+")
	{
		return (substr $chromosome_sequence{$gff_feature{$id}->{Chr}}->{Seq}, $gff_feature{$id}->{Start}-1, $gff_feature{$id}->{End}-$gff_feature{$id}->{Start}+1);
	}
	else
	{
		my $sequence = substr $chromosome_sequence{$gff_feature{$id}->{Chr}}->{Seq}, $gff_feature{$id}->{Start}-1, $gff_feature{$id}->{End}-$gff_feature{$id}->{Start}+1;
		#		my $reverse_complementary = reverse $sequence; 
		#		$reverse_complementary =~ tr/[A, T, C, G]/[T, A, G, C]/;
		my $reverse_complementary = reverse $sequence;
		$reverse_complementary =~ tr/[A, T, C, G, a, t, c, g]/[T, A, G, C, t, a, g, c]/;
		return $reverse_complementary;
	}
}

sub get_promoter_sequence {
	# body...
	my $id = shift; 
	
	if ($gff_feature{$id}->{Strand} eq "+")
	#on the watson strand, looking upstream for promoter
	{
		my $promoter_start_site = $gff_feature{$id}->{Start}-500;
		#Now look for the last gene
		
		my $last_gene = search_last_gene($id);
		
		#foreach my $key (sort {$a cmp $b} keys %gff_feature) {
		#	if ($gff_feature{$key}->{Gene_order} eq $last_gene_order)
		#	{
		#		$last_gene = $key;
		#	}
		#}
		#print $last_gene."\n";
		
		#Only consider the last gene boundary when both genes are located in the same chromosome
		if ($last_gene && ($gff_feature{$id}->{Chr} eq $gff_feature{$last_gene}->{Chr})) {
			my $last_gene_boundary = $gff_feature{$last_gene}->{End};
			if ($last_gene_boundary > $promoter_start_site) {
				$promoter_start_site = $last_gene_boundary;
			}
		}
		
		#print $promoter_start_site."\n";
		return substr($chromosome_sequence{$gff_feature{$id}->{Chr}}->{Seq}, $promoter_start_site-1, $gff_feature{$id}->{Start}-$promoter_start_site); 
	}
	else
	{
		#on the crick strand, looking downstream for promoter
		#remember to reverse complement the promoter sequence
		my $promoter_start_site = $gff_feature{$id}->{End}+500;
		
		#now look for the next gene
		#my $next_gene_order = $gff_feature{$id}->{Gene_order}+1;
		
		my $next_gene = search_next_gene($id);
		
		#foreach my $key (sort {$a cmp $b} keys %gff_feature) {
		#	if ($gff_feature{$key}->{Gene_order} eq $next_gene_order)
		#	{
		#		$next_gene = $key;
		#	}
		#}
		#print $next_gene."\n";
		#print $gff_feature{$next_gene}->{Start}."\n";
		#print $gff_feature{$id}->{End}."\n";
		#print $promoter_start_site."\n";
		#only consider the next gene boundary when both genes are in the same chromosome
		if ($gff_feature{$id}->{Chr} eq $gff_feature{$next_gene}->{Chr}) {
			# body...
			my $next_gene_boundary = $gff_feature{$next_gene}->{Start}; 
			if ($next_gene_boundary < $promoter_start_site) {
				$promoter_start_site = $next_gene_boundary;
			}
		}
		my $promoter_sequence = substr($chromosome_sequence{$gff_feature{$id}->{Chr}}->{Seq}, $gff_feature{$id}->{End}, $promoter_start_site-$gff_feature{$id}->{End});
		$promoter_sequence = reverse $promoter_sequence;
		$promoter_sequence =~ tr/[A, T, C, G, a, t, c, g]/[T, A, G, C, t, a, g, c]/;
		return $promoter_sequence;
	}
	
}

sub get_3UTR_sequence {
	# body...
	my $id = shift; 
	if ($gff_feature{$id}->{Strand} eq "+")
	#on the watson strand, looking downstream for 3'UTR
	{
		my $UTR_end_site = $gff_feature{$id}->{End}+200;
		#Now look for the next gene
		#my $next_gene_order = $gff_feature{$id}->{Gene_order}+1;
		
		my $next_gene = search_next_gene($id);
		
		#foreach my $key (sort {$a cmp $b} keys %gff_feature) {
		#	if ($gff_feature{$key}->{Gene_order} eq $next_gene_order)
		#	{
		#		$next_gene = $key;
		#	}
		#}
		#print $next_gene."\n";
		
		#Only consider the next gene boundary when both genes are located in the same chromosome
		if ($next_gene && ($gff_feature{$id}->{Chr} eq $gff_feature{$next_gene}->{Chr})) {
			my $next_gene_boundary = $gff_feature{$next_gene}->{Start};
			if ($next_gene_boundary < $UTR_end_site) {
				$UTR_end_site = $next_gene_boundary;
			}
		}
		
		#print $promoter_start_site."\n";
		return substr($chromosome_sequence{$gff_feature{$id}->{Chr}}->{Seq}, $gff_feature{$id}->{End}, $UTR_end_site-$gff_feature{$id}->{End});
	}
	else
	{
		#on the crick strand, looking upstream for 3'UTR
		#remember to reverse complement the promoter sequence
		my $UTR_end_site = $gff_feature{$id}->{Start}-200;
		
		#now look for the last gene
		#my $last_gene_order = $gff_feature{$id}->{Gene_order}-1;
		
		my $last_gene = search_last_gene($id);
		
		#foreach my $key (sort {$a cmp $b} keys %gff_feature) {
		#	if ($gff_feature{$key}->{Gene_order} eq $last_gene_order)
		#	{
		#		$last_gene = $key;
		#	}
		#}
		#print $last_gene."\n";
		#print $gff_feature{$last_gene}->{End}."\n";
		#print $gff_feature{$id}->{Start}."\n";
		#print $UTR_end_site."\n";
		#only consider the last gene boundary when both genes are in the same chromosome
		if ($gff_feature{$id}->{Chr} eq $gff_feature{$last_gene}->{Chr}) {
			# body...
			my $last_gene_boundary = $gff_feature{$last_gene}->{End}; 
			if ($last_gene_boundary > $UTR_end_site) {
				$UTR_end_site = $last_gene_boundary;
			}
		}
		my $UTR_sequence = substr($chromosome_sequence{$gff_feature{$id}->{Chr}}->{Seq}, $UTR_end_site, $gff_feature{$id}->{Start}-$UTR_end_site);
		$UTR_sequence = reverse $UTR_sequence;
		$UTR_sequence =~ tr/[A, T, C, G, a, t, c, g]/[T, A, G, C, t, a, g, c]/;
		return $UTR_sequence;
	}	
}

sub check_restriction_sites {
	#BsaI: GGTCTC
	#BsmBI: CGTCTC
	my $sequence_to_check = shift; 
	my $restriction_status; 
	if ($sequence_to_check =~ m/ggtctc|gagacc/i) {
		$restriction_status .= "BsaI, ";
	}
	if ($sequence_to_check =~ m/cgtctc|gagacg/i) {
		$restriction_status .= "BsmBI";
	}
	if (defined($restriction_status)) {
		return $restriction_status; 
	}
	else {
		$restriction_status = "no BsaI or BsmBI site";
		return $restriction_status;
	}
}


sub tm_cal {
	my $sequence = shift; 
	my $gc_counts = gc_content_cal($sequence);
	my $tm; 
	if (length($sequence)<=14)
	{
		$tm = 4*length($sequence)*$gc_counts + 2*length($sequence)*(1-$gc_counts);
		
	}
	else
	{
		$tm = 64.9 + 41*($gc_counts*length($sequence)-16.4)/length($sequence);
	}
	return $tm = int($tm);
}


sub gc_content_cal {
	my $sequence = shift; 
	my $gc_counts = ($sequence =~ tr/G|g|C|c//)/length($sequence);
	return $gc_counts;
}

sub reverse_complement {
	my $sequence = shift; 
	$sequence = reverse $sequence;
	$sequence =~ tr/[A, T, C, G, a, t, c, g]/[T, A, G, C, t, a, g, c]/;
	return $sequence;
}



sub design_golden_gate_primers {
	my $id = shift;
	my %primer_sequences = (); 
	print "\nDESIGNING PRIMERS FOR GENE $id\n";
	print "=================================\n";
	my $promoter_5primer = "agcgtgGGTCTCgGGCT"; 
	my $promoter_3primer = "GATGtGAGACCcagcac";
	my $ORF_5primer = "agcgtgGGTCTCaG"; 
	my $ORF_3primer = "TAGCcGAGACCcagcac";
	my $UTR_5primer = "agcgtgGGTCTCtTAGC";
	my $UTR_3primer = "CCTCcGAGACCcagcac";  #length is 13, so the priming region should not be more than 47


	my $promoter_sequence = get_promoter_sequence($id);
	my $ORF_sequence = get_orf_sequence($id);
	my $UTR_sequence = get_3UTR_sequence($id);

	

#	my @introns = check_introns($id); 

	$primer_sequences{$id}->{Promoter} = get_promoter_sequence($id);
	$primer_sequences{$id}->{ORF} = get_orf_sequence($id);
	$primer_sequences{$id}->{UTR} = get_3UTR_sequence($id);




		#take initial 18bp and extend it up to 25bp, choose the one which has TM closet to 60℃. 
		#add 5' leading sequence or 3' following sequence to F and R primers respetively
		#make sure the R primer is reverse complementary 
		#output Ape file (optional)
		my $diff = 99;
		my $current_best_pos = 18; 
		for (my $i = 18; $i<=30; $i++)
		{
			#print substr($promoter_sequence, 0, $i)." \n";
			my $promoter_sub_seq = substr($promoter_sequence, 0, $i);
			my $tm_promoter_sub_seq = tm_cal($promoter_sub_seq);
			#print $promoter_sub_seq."  ".length($promoter_sub_seq)." ".$tm_promoter_sub_seq." ".abs($tm_promoter_sub_seq-60)."\n";
			if (abs($tm_promoter_sub_seq-60) < $diff)
			{
				$diff = abs($tm_promoter_sub_seq - 60);
				$current_best_pos = $i;
			}
		}
		my $promoter_F_primer = $promoter_5primer.substr($promoter_sequence, 0, $current_best_pos);
		$primer_sequences{$id}->{promoter_F_primer} = $promoter_F_primer; 
		print $id."_P_F\t".$promoter_F_primer."\tTM\t".tm_cal(substr($promoter_sequence, 0, $current_best_pos))."\n";
		my $primer_name = $id."_P_F";
		
		$forward_primers{$primer_name} = $promoter_F_primer;
		
		#design $promoter_R_primer
		$diff = 99;
		$current_best_pos = 18; 
		
		for (my $i = 18; $i<=30; $i++)
		{
			#print substr($promoter_sequence, -$i)." ".length(substr($promoter_sequence, -$i))." \n";
			my $promoter_sub_seq = substr($promoter_sequence, -$i);
			my $tm_promoter_sub_seq = tm_cal($promoter_sub_seq);
			#print $promoter_sub_seq." ".$tm_promoter_sub_seq."\n";
			if (abs($tm_promoter_sub_seq - 60) < $diff)
			{
				$diff = abs($tm_promoter_sub_seq - 60);
				$current_best_pos = $i;
			}
		}
		my $promoter_R_primer = substr($promoter_sequence, -$current_best_pos).$promoter_3primer;
		$promoter_R_primer = reverse_complement($promoter_R_primer);
		$primer_sequences{$id}->{promoter_R_primer} = $promoter_R_primer; 
		
		print $id."_P_R\t".$promoter_R_primer."\tTM\t".tm_cal(substr($promoter_sequence, -$current_best_pos))."\n";		
		$primer_name = $id."_P_R";
		$reverse_primers{$primer_name} = $promoter_R_primer; 
		
		#design $orf_F_primer
		$diff = 99; 
		$current_best_pos = 18; 
		
		for (my $i = 18; $i <= 30; $i++)
		{
			my $orf_sub_seq = substr($ORF_sequence, 0, $i);
			my $tm_orf_sub_seq = tm_cal($orf_sub_seq);
			#print $orf_sub_seq." ".length($orf_sub_seq)." ".$tm_orf_sub_seq."\n";
			if (abs($tm_orf_sub_seq - 60) < $diff)
			{
				$diff = abs($tm_orf_sub_seq -60);
				$current_best_pos = $i;			
			}
		}
		my $orf_F_primer = $ORF_5primer.substr($ORF_sequence, 0, $current_best_pos);
		$primer_sequences{$id}->{orf_F_primer} = $orf_F_primer; 
		print $id."_G_F\t".$orf_F_primer."\tTM\t".tm_cal(substr($ORF_sequence, 0, $current_best_pos))."\n";
		$primer_name = $id."_G_F";
		$forward_primers{$primer_name} = $orf_F_primer;
		
		#design $orf_R_primer
		$diff = 99; 
		$current_best_pos = 18; 
		
		for (my $i = 18; $i <= 30; $i++)
		{
			my $orf_sub_seq = substr($ORF_sequence, -$i);
			#print $orf_sub_seq." ".length($orf_sub_seq)."\n";
			my $tm_orf_sub_seq = tm_cal($orf_sub_seq);
			if (abs($tm_orf_sub_seq - 60) < $diff)
			{
				$diff = abs($tm_orf_sub_seq);
				$current_best_pos = $i;
			}			
		}
		my $orf_R_primer = substr($ORF_sequence, -$current_best_pos).$ORF_3primer;
		$orf_R_primer = reverse_complement($orf_R_primer);
		$primer_sequences{$id}->{orf_R_primer} = $orf_R_primer; 
		print $id."_G_R\t".$orf_R_primer."\tTM\t".tm_cal(substr($ORF_sequence, -$current_best_pos))."\n";
		$primer_name = $id."_G_R";
		$reverse_primers{$primer_name} = $orf_R_primer;
		
		#design $utr_F_primer
		$diff = 99;
		$current_best_pos = 18; 
		
		for (my $i = 18; $i <= 30; $i++)
		{
			my $utr_sub_seq = substr($UTR_sequence, 0 , $i); 
			my $tm_utr_sub_seq = tm_cal($utr_sub_seq);
			if (abs($tm_utr_sub_seq - 60) < $diff)
			{
				$diff = abs($tm_utr_sub_seq - 60);
				$current_best_pos = $i;
			}
		}
		my $utr_F_primer = $UTR_5primer.substr($UTR_sequence, 0, $current_best_pos);
		$primer_sequences{$id}->{utr_F_primer} = $utr_F_primer; 
		print $id."_T_F\t".$utr_F_primer."\tTM\t".tm_cal(substr($UTR_sequence, 0, $current_best_pos))."\n";
		$primer_name = $id."_T_F";
		$forward_primers{$primer_name} = $utr_F_primer;
		
		#design $utr_R_primer
		$diff = 99;
		$current_best_pos = 18;
		
		for (my $i = 18; $i <= 30; $i++)
		{
			my $utr_sub_seq = substr($UTR_sequence, -$i);
			my $tm_utr_sub_seq = tm_cal($utr_sub_seq);
			if (abs($tm_utr_sub_seq - 60) < 60)
			{
				$diff = abs($tm_utr_sub_seq - 60);
				$current_best_pos = $i;
			}
		}
		my $utr_R_primer = substr($UTR_sequence, -$current_best_pos).$UTR_3primer;
		$utr_R_primer = reverse_complement($utr_R_primer);
		$primer_sequences{$id}->{utr_R_primer} = $utr_R_primer; 
		print $id."_T_R\t".$utr_R_primer."\tTM\t".tm_cal(substr($UTR_sequence, -$current_best_pos))."\n";
		$primer_name = $id."_T_R";
		$reverse_primers{$primer_name} = $utr_R_primer;	
		
		return %primer_sequences;			
}




