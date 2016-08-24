use strict;
use warnings;



##list blast databases
my $wd = '/Users/farman/blast_databases/';
my @files = <$wd*.fasta>;

foreach my $file (@files) {   ##for each blast database	

system("blastn -query primer3_primers2-17-16.fasta -outfmt '6 qseqid sseqid qlen length pident qstart qend sstart send' -db $file -word_size 16 >>outnewprimers.txt")

#print "blastn -query test_primers.fasta -perc_identity 100 -db $file";


}


