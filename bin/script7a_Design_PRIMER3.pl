##!/usr/bin/env perl
##Design primer3 file based on FASTA files.  execute in directory with final FASTA files.
#Bradford Condon, University of Kentucky

use strict;
use warnings;
use Bio::SeqIO;

my @files = <*.fasta>;

 foreach my $input_file (@files) { 
	my $seqio = Bio::SeqIO->new(-file=>$input_file, '-format' => 'Fasta');
	my $output_file = "$input_file.primer3.txt";
	open( OUT, ">$output_file" )
	or die "error: cannot open $output_file for writing: $!\n";

	while( my $seqObject = $seqio->next_seq ){
	
	# get sequence information
	my $id  = $seqObject->id;
	my $seq = $seqObject->seq;

	print OUT "PRIMER_SEQUENCE_ID=$id\n";
	print OUT "SEQUENCE_TEMPLATE=$seq\n";
	print OUT "PRIMER_PRODUCT_SIZE_RANGE=370-400\n";
	print OUT "PRIMER_MAX_SIZE=23\n";
	print OUT "PRIMER_NUM_RETURN=1\n";
	print OUT "=\n";



	}
}