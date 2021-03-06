#blasts a single FASTA file against a directory of blast databases.
#output is piped to a single file (.blast extension)
use strict;
use warnings;
use Getopt::Long;

my $usage = "$0 -i Fasta to BLAST  \n"; 

my $input_file;

#read user options
GetOptions(
	"i=s"	=> \$input_file,
);

##list blast databases
my $wd = '/Users/chet/uky/blast_dbs/'; #/home/bradford.condon/DiagLockPick/blast_dbs/

my @files = <$wd*.fasta>;

foreach my $file (@files) {   ##for each blast database	

system("blastn -query $input_file -db $file -max_target_seqs 2   >> $input_file.blast");

}


