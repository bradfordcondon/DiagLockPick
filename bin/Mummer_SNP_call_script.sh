###
#January 2016
#Bradford Condon
#Farman Lab
#University of Kentucky
###Goal of script:
#Automate the calling of SNPs by MUMmer v.3.2.0

#mkdir header_date

mkdir MUMmerSNPs

for file in *.fasta
do
###
#mum option flag specifies use only 1:1 matches for initial matches
nucmer --mum WBKY11_scaffolds.fasta $file 

#delta file is out.delta
#run dnadiff to call SNPs
dnadiff -d out.delta
##
#Note that in the output, the number of repeat alignments which cover this reference position [Q] number of 
#repeat alignments which cover this query position 
#to exclude these, we'll use the C option

show-snps -C out.delta >  out_SNPs.txt

mv out_SNPs.txt MUMmerSNPs/$file.msnps.txt

done