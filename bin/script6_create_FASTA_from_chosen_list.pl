#!/usr/bin/env perl
#simple Bioperl script to extract seqs from FASTA file within a certain range.
# Program: order scaffolds

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Data::Dumper;
# usage
my $usage = "$0 -i|input fasta_input_files -l|list list of scaffold and ranges \n";

# global values
my $input_file;
my $list_file;
# read user options
GetOptions(
	"i|input=s" => \$input_file,
	"l|list=s" => \$list_file,

);

# check for user parameters
if( !$input_file ){
	die $usage;
}
if( !$list_file ){
	die $usage;
}

#open list file

my %listhash;

open(my $fh, "<", $list_file)
			or die "couldnt open '$list_file' $!";
						
			while (<$fh>){
				chomp;
				my @split = split(/[_-]/);


				$listhash{$split[0]}{$split[1]} =$split[2];
				
					
			}

print Dumper(\%listhash);




# open input fasta file
unless( -e $input_file ){
	die "error: cannot find fasta input file $input_file\n";
}
my $input = Bio::SeqIO->new (-file => "<$input_file", '-format' => 'Fasta')
	or die "error: failure opening fasta $input_file for reading: $!\n";

# make temporary fasta file and print header line
my $output_file = "remaining_loci.fasta";
open( OUT, ">$output_file" )
	or die "error: cannot open $output_file for writing: $!\n";
# step through sequences in input fasta file
while( my $seqObject = $input->next_seq ){
	
	# get sequence information
	my $id  = $seqObject->id;
	$id =~ s/scaffold/s/;

	my $seq = $seqObject->seq;

		foreach my $k (keys %listhash) {  

	
				if ($k eq $id) {  

	my @seqarray = split('', $seq);
			
#for each start, print out the range fmor start to end
			foreach my $start (keys $listhash{$k}) {
		

					my $end = $listhash{$k}{$start};

					print "extracting $k from $start to $end\n";
				
					my @slice = @seqarray[$start..$end];

					my $outseq = join('', @slice);

					  print OUT ">$id","_","$start-$end\n", "$outseq\n";  

					  }
				}
		}
					
	}

