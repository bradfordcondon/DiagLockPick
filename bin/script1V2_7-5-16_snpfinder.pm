 #! /usr/bin/perl -w
	#this script goes through SNP output.  It counts the number of times a scaffold has SNPs compared to another genome.
	#It outputs the number of genomes that have SNPs on that scaffold. 
	#updated July 29 2016 with no clade distinction, but working on MUMMer.
	
	use strict;
	use warnings;
	use Data::Dumper;  #for debugging
	use List::Util qw(max);  #for getting maximum hash value for sliding window
	use List::Util qw(sum);
	use Getopt::Long;

	my $usage = "$0 -d|directory directory -p|conserved TRUE/FALSE generate conserved region file \n 
	-c|cladelist 3 column clade list -r TRUE/FALSE generate report file -i mummer (default = FALSE) \n";
	my $directory;
	my $cladelist;
	my $reportprint = "FALSE";
	my $conservedprint = "FALSE";
	my $windowsize;
	my $stepsizemaster;
	my $minsnpsite;
	my $mummer;
	my $maxPrimerSnp;
	
	# read user options
	GetOptions(
		"d|directory=s" => \$directory,
		"c|cladelist=s" => \$cladelist,
		"r|report=s" =>\$reportprint,
		"p|conserved" =>\$conservedprint,
		"w|windowsize=s" =>\$windowsize,
		"s|stepsize=s" =>\$stepsizemaster,
		"m|minsnpsite=s" =>\$minsnpsite,
		"i|mummer=s" =>\$mummer,
		"x|maxPrimerSNP=s" =>\$maxPrimerSnp
			);

		# check for user parameters	
			if( !$directory ){
				$directory ='/Users/chet/uky/blast_dbs/MUMmerSNPs/'; print "no directory specified (-d), using /Users/chet/uky/blast_dbs/MUMmerSNPs\n";
			}
		 if(!$cladelist ) {
	 $cladelist = '/Users/chet/uky/list_mummers.txt'; print "no cladelist specified (-c), using default list.txt\n";
	 		}
 if(!$windowsize ) {
	 $windowsize = '400'; print "no windowsize specified(-w), using default 400\n";
	 		}
if(!$stepsizemaster ) {
	 $stepsizemaster = '50'; print "no step size specified(-w), using default 50\n";
	 		}		
if (!$minsnpsite) {
	$minsnpsite = '0'; print "no min SNP site/window set (-m), using 0 default\n";
}
if (!$mummer) {
	$mummer = 'TRUE'; print "no MUMmer option (-i) set, using TRUE \n";
}
if (!$maxPrimerSnp) {
	$maxPrimerSnp = 40; print "no max primer SNP -x set, using default = 40\n";
}
	my %masterhash;  #declare global variables. 
	# Master hash is the main hash, exclusion hash is the SNPs to exlclude
	my %exclhash;
	my %allhash;
	my %allnucs;

	my @files = <$directory*>;  #separate reference files
	print $files[0];
		
 ###################################################                	   
	### Load a list of species, host, clade. 
	### Assign them to in or out clade arrays based on if we want them to group with or outside of the reference
 ###################################################                	   

my %cladetracker;
	
open(my $fh, "<", $cladelist)
or die "couldnt open '$cladelist' $!";		
	while (<$fh>){
	chomp;
	my @split = split(/\s+/);
	my $name = $split[0];
	my $clade = $split[2];
	$name = $directory.$name;
	
	$cladetracker{$name} = $clade;					
}
close $fh;


###################################################                	   
	### read in masks as hash
 ###################################################    

my %maskhash;
#my $mask_file = "WBKYmaplist.txt";
#print "reading mask file\n";
#%maskhash = build_SNP_NUC_hash($mask_file, \%maskhash);
#close $mask_file;
#print Dumper(\%maskhash);  #currently empty

 ###################################################                	   
	####	Go through each SNP report and build a hash based on 
 ###################################################                	   
			
	foreach my $file (@files) {   ##open each SNP file
	print "working on $file \n!";
			my $element = $file;
######################
# Build hash of each nucleotide change at the SNP
######################			
		%allnucs = build_SNP_NUC_hash($file, \%allnucs);
	   	%allhash = build_SNP_hash($file, \%allhash);	 	 				 		
			}								  			 		 
 ###################################################                	   
 ###################################################                	   


###################################################                	   
	### Load in list of scaffold lengths.
 ###################################################    


open(my $length_file, "<", "scaffold_summaries.txt")
or die "couldn't open length data file\n";

my %lengthinfo;
while (<$length_file>){
	chomp;
	my @split = split(/\s+/);
	my $key = $split[0];
	my $value = $split[1];
	$lengthinfo{$key} = $value;
}

print Dumper(\%lengthinfo);


######################
#count unique character instances in each key
######################

my %reporthash;	

#k is scaffold, k2 is nucleotide position.  

 	 foreach my $k (keys %allnucs) {   
	   		foreach my $k2 (keys $allnucs{$k}) {	

					   		my $k3= $allnucs{$k}{$k2};	#k3 is a mash of nucleotides 
					   		
					   		 my $count =0;
     						 my %hashofuniq; 
					   		
					   			 my @nucstring = split ('',$k3);  #splits scalar of all nucleotides at a site.
					   			 @hashofuniq{@nucstring} = 1;   #essentially counts # unique nucleotides
					   			 my $length =  scalar @nucstring;   #returns a length scalar for the number of elements in above array. 
      foreach  (keys %hashofuniq){
   					     $count++; 
     							 }

     							 $reporthash{$k}{$k2}{"Numnucs"}=$count;
     							 $reporthash{$k}{$k2}{"NumSNPs"}=$length;
												}
									}
	###################################################                	   
 ###################################################                	   
		

#####
#Output  nucleotide Report from report hash

if ($reportprint eq "TRUE" ) {   ###only print report if desired
		print "outputting report file because -r =TRUE\n";

	my   $reportfile = 'nucleotide_diversity_report.txt';
	open (my $reportout, '>', $reportfile) or die "could not open $reportfile $!";
	print $reportout "Scaffold\tTotal SNP sites\tAvg SNPs per base\t Average nucleotides per SNP\n";

	foreach my $k (sort keys %reporthash) {   

		my @tracknumSNPs;
		my @tracknumnucs;

		foreach my $k2 (keys $reporthash{$k}) {	
	   			my $Numnuc = $reporthash{$k}{$k2}{"Numnucs"};
	   			my $numSNP = $reporthash{$k}{$k2}{"NumSNPs"};
	   			push @tracknumnucs, $Numnuc;
	   			push @tracknumSNPs, $numSNP;
		}

   		my $avgnucs =  sum(@tracknumnucs)/@tracknumnucs;
   		my $avgsNPS =  sum(@tracknumSNPs)/@tracknumSNPs;
   		my $count =  @tracknumnucs;

   		print $reportout "$k\t$count\t$avgsNPS\t$avgnucs\n";
	}		
		}	


 ###################################################                	   
	###Build file of conserved regions.
 ###################################################                	   


if ($conservedprint eq "TRUE") {   ###only print conserved file if desired
			print "outputting conserved region file because -p =TRUE\n";

my   $conserved = 'conserved.txt';
	    open (my $consfile, '>', $conserved) or die "could not open $conserved $!";
	
	
	my %NoSNP= sliding_window(\%allhash, \%allhash, "10", "30");


	#output: a hash of SNP counts for 30bp windows.
	
	#for conserved regions, I want to step through allhash, and check against the allhash
  
	  foreach my $k (sort keys %NoSNP) {
  	foreach my $k2(sort {$a<=>$b} keys $NoSNP{$k}) {			
  			my $k3= $NoSNP{$k}{$k2}{"count"};
  			my $end = $NoSNP{$k}{$k2}{"end"};
  
  	if ( $k3 == 0) {  ##only show ranges with a zero
	  print $consfile "$k\t$k3\t$end\n" ; ###print out keys at end

					}
													}
									}
		
	 close $consfile;

}
	 
	     ######Remove exclusion hash from allhash , creating allhash
		#######
		
	
 ###At this point we may want to only show SNPs with more than X counts- ie those SNPs showing up in more than X of the query genomes.
###As a cautious starting point, lets exclude SNPs that show up in <50% of the other genomes
	   

	   #for now: minimum number of other genomes = 1.  Not filtering, essentially.

	
my $outfile = 'out.txt';
open (my $oh, '>', $outfile) or die "could not open $outfile $!";
	
  
		
	print "performing final window analysis\n";

#step size then windowsize
	my %SNP_rich_windows= Diversity_sliding_window(\%allhash, \%allhash, $stepsizemaster, $windowsize);	  
	      	      	      	     


####



	      	      	print $oh "Scaffold\tstart\tend\tNumsnps\tdivscore\tavg_participants\tmaskcount\tprimercount\n";

	foreach my $k (sort keys %SNP_rich_windows) {
		print "working on $k\n";

  		foreach my $k2( sort keys $SNP_rich_windows{$k}) { 

  					my $snpcount =$SNP_rich_windows{$k}{$k2}{"count"};
  					my $end = 	$SNP_rich_windows{$k}{$k2}{"end"};
  					my $predivscore =$SNP_rich_windows{$k}{$k2}{"diversity"};
 					my $maskfinal = $SNP_rich_windows{$k}{$k2}{"mask"};
 					my $participants = 	$SNP_rich_windows{$k}{$k2}{"participants"};
 					my $primerCount = $SNP_rich_windows{$k}{$k2}{"primerCount"};	

 					my $finalparts =0;
 					my $divscore = 0;
 					unless ($snpcount == 0 ){
 					 $divscore = $predivscore/$snpcount;
 					  $finalparts = $participants/$snpcount;

 				}

 					if ($snpcount >= $minsnpsite ){
 						if ($maskfinal <= ($windowsize/2) ){

 							if ($primerCount <= $maxPrimerSnp) {
		
					  	print $oh "$k\t$k2\t$end\t$snpcount\t$divscore\t$finalparts\t$maskfinal\t$primerCount\n"; ###print out keys at end:	  	
					}
					}
				}
		}
	}			
close $oh;
			
			
	###Subroutines######
			
			
	########
	########
	#Sliding window subroutine.  
	#This Subroutine takes two 2-layer hash.  It sorts them by scaffolds, and builds an array of nucleotides for each scaffold.
	#It  then slides along at a defined stepsize and window size.  It scrolls through the first (%hash) and checks against teh second %hashcheck
	#If there is a value in the input hash found as it slides along, it counts it.
	#It will output a 2 level hash: scaffold, range,  and count.
			
			##Input is hash reference, step size, windowsize.
			
			
				sub sliding_window {  
				
   				 my %hash = %{shift()};
   				 my %hashcheck = %{shift()};
   				 my $stepsize = shift;
   				 my $windowsize = shift; 
				 my %final;
				
	foreach my $k (sort { $hash{$a} <=> $hash{$b} } keys %hash){   ##sort the input hash																					
	
					my @locsort;	   #define the sorted array for each scaffold
					foreach my $k2 (keys ($hash{$k})) {	     
						push (my @loc, $k2);  #build array of nucleotides
						@locsort = sort {$b<=>$a} @loc;   #sort the array of nucleotides from biggest to smallest.  Backwards because we want to know what the max value is below.
														}
	
				my $p = 1;  #p is my scaffold bookmark.  We will always start at 1 for each scaffold
				my $max = $lengthinfo{$k};  #end of the scaffold. 
	
				for ($p; $p< $max; $p= $p+$stepsize) {   #this moves our window
	
						my $count = 0;				#we're counting SNPs- restart counter at 0 for each window
						my $i = $p;					#start window at scaffold bookmark	
						my $winsize = $windowsize+$i;			#our window range is our start to start + windowsize
	
					for ($i; $i <= $winsize; $i++) {	#check each nucleotide, 1 by 1, in our range, for a hash key		

			$count++	if exists $hashcheck{$k}{$i};  #if we found a hash key, count it.  
					}
								
			my $end = $i;
			if ($end > $max){
				$end = $max;
			}  
			my $start = $i-$windowsize;
			  
			$final{$k}{$start}{"count"} = $count;  #assing range to count value	
			$final{$k}{$start}{"end"}	= $end;							
													}  #done looping through  windows	
	  }	#done looping through scaffolds	
				return (%final);				
				}
				


#######			
#######Sub filter window.  Improves upon above sliding window script.  
#######

				sub Diversity_sliding_window {  
				
   				 my %hash = %{shift()};
   				 my %hashcheck = %{shift()};
   				 my $stepsize = shift;
   				 my $windowsize = shift; 
				 my %final;
				
	foreach my $k (sort { $hash{$a} <=> $hash{$b} } keys %hash){   ##sort the input hash																					
	
					my @locsort;	   #define the sorted array for each scaffold
					foreach my $k2 (keys ($hash{$k})) {	     
						push (my @loc, $k2);  #build array of nucleotides
						@locsort = sort {$b<=>$a} @loc;   #sort the array of nucleotides from biggest to smallest.  Backwards because we want to know what the max value is below.
														}
	
				my $p = 1;  #p is my scaffold bookmark.  We will always start at 1 for each scaffold
				my $max = $lengthinfo{$k};  #end of the scaffold. 
	
				for ($p; $p< $max; $p= $p+$stepsize) {   #this moves our window
	
						my $count = 0;	#we're counting SNPs- restart counter at 0 for each window
						my $primerCount = 0;
						my $windivscore = 0;   #set the window's diversity score to 0.
						my $i = $p;					#start window at scaffold bookmark	
						my $winsize = $windowsize+$i;			#our window range is our start to start + windowsize
						my $maskcount = 0;
						my $partscore = 0;
					for ($i; $i <= $winsize; $i++) {	#check each nucleotide, 1 by 1, in our range, for a hash key		

			$count++	if exists $hashcheck{$k}{$i};  #if we found a hash key in the hash we are counting things, count it.  
						if ($winsize - $i <= 20  || $winsize - $i > $windowsize - 20){#Checking here for exact conserved start/end for loci.
														#this will ensure that primers are likely to amplify.
											$primerCount++	if exists $hashcheck{$k}{$i}; 			
						}

						if (exists $reporthash{$k}{$i}) {

					my $toadd =  $reporthash{$k}{$i}{"Numnucs"};  #score is simply the number of unique nucleotides- 1-3 (A, T, G, C- reference is assumed as 4th nucleotdie).
					my $participants =  $reporthash{$k}{$i}{"NumSNPs"};  #Number of strains participating in each site
					$windivscore = $windivscore+$toadd;
					$partscore = $partscore+$participants;
						}				
						if (exists $maskhash{$k}{$i}){
							$maskcount++;
						}					
					}
								
			my $end = $i;  #build keystring as name range
			my $start = $i-$windowsize;

			#Don't want windows to be longer than the scaffold- fix that here.
			if ($end > $max){
				$end = $max;
			}  
						
			$final{$k}{$start}{"primerCount"} = $primerCount;		
			$final{$k}{$start}{"participants"} = $partscore;
			$final{$k}{$start}{"count"} = $count;  #assign range to count value	
			$final{$k}{$start}{"end"}	= $end;		
			$final{$k}{$start}{"diversity"}	= $windivscore;			#assign diversity score
			$final{$k}{$start}{"mask"} = $maskcount;
													}  #done looping through  windows	
	  }	#done looping through scaffolds	
				return (%final);				
				}
				


####OPEN SNP FILE AND BUILD HASH BASED ON COUNTS

	sub build_SNP_hash { 
	
				 my $file = shift;
				 my %outhash = %{shift()};
				my $linecount = 0;
	open (my $fh, "<", $file)
	  	  or die "couldn't open '$file' $!"; 		

	  	  if ($mummer eq 'TRUE' ) {
			 while (<$fh>) {
			 if ($linecount > 5) {#remove first 5 lines
				 chomp;
				 $linecount++;
				 }
				else {
	        	chomp ;
	        	next unless length;	
			    my @split = split(/\s+/); 
	        	next unless $split[11];	  	  	
			    my $queryloc= $split[11];
		  	    my $hitname= $file;
		  	    my $quernuc = $split[0];
		  	    my $querbase = $split[1];
		  	    my $hitbase = $split[2];

			  	 unless ($hitbase eq 'N' ||$hitbase eq '.' || $querbase eq 'N' || $querbase eq '.'){
		  	     	++$outhash{$queryloc}{$quernuc};  ##keep track of the number of SNPs at each basepair, for each scaffold
		  	    	}
		  	    	}
		  	}
			return (%outhash); 	
			}
				
	  	  else {
	   		 while (<$fh>) {
			    chomp ;
			    next unless length;
			    my @split = split(/\s+/);   
			    my $queryloc = $split[0];
			    my $hitname = $split[1];
			    my $quernuc = $split[2];
				++$outhash{$queryloc}{$quernuc};  ##keep track of the number of SNPs at each basepair, for each scaffold
			      } 						
	   		return (%outhash);
	   		 } 
 }


####OPEN SNP FILE AND BUILD HASH OF NUCLEOTIDE OCCURRENCES



 
sub build_SNP_NUC_hash {

		 my $file = shift;
		 my %outhash = %{shift()};
			open (my $fh, "<", $file)  or die "couldn't open '$file' $!"; 
			my $linecount=1;	
				#separate parsing for mummer or non-mummer
			if ($mummer eq 'TRUE' ) {
				 while (<$fh>) {
				 if ($linecount > 5) {#remove first 5 lines
				 chomp;
				 $linecount++;
				 }
				 else{		 
		        	chomp ;
		        	next unless length;	
				    my @split = split(/\s+/); 
		        	next unless $split[11];	  	  	
				    my $queryloc= $split[11];
			  	    my $hitname= $file;
			  	    my $quernuc = $split[0];
			  	    my $hitbase = $split[2];
			  	    my $querbase = $split[1];
			  	    unless ($hitbase eq 'N' ||$hitbase eq '.' || $querbase eq 'N' || $querbase eq '.'){
	  	    		 $outhash{$queryloc}{$quernuc} .="$hitbase";  ##keep track of the number of SNPs at each basepair, for each scaffold
	  	   			}
	  	   		}
	  	    	}
				return (%outhash);
  	 		}	
		 else {   					
			while (<$fh>) {
		    chomp ;
		    next unless length;
		    my @split = split(/\s+/);   
		    my $queryloc = $split[0];
		    my $hitname = $split[1];
		    my $quernuc = $split[2];
		    my $hitnuc = $split[3];
			my $hitbase = $split[5];
			$outhash{$queryloc}{$quernuc} .= "$hitbase";
				}
	  	return (%outhash);
		}
}
				
				