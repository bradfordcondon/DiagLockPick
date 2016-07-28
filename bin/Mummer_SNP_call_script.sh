###
#January 2016
#Bradford Condon
#Farman Lab
#University of Kentucky

###
###Goal of script:
#Automate the calling of SNPs by MUMmer v.3.2.0

###Step one: navigate to and set up working directory
cd /Users/chet/uky/mummerSNPs_Jan2016
#mkdir header_date

###
#mum option flag specifies use only 1:1 matches for initial matches
nucmer --mum Arcadia_merged.fasta Br80_merged.fasta 

#delta file is out.delta
#run dnadiff to call SNPs
dnadiff -d out.delta
##
#Note that in the output, the number of repeat alignments which cover this reference position [Q] number of 
#repeat alignments which cover this query position 
#to exclude these, we'll use the C option

show-snps -C out.delta >try1snps.txt

