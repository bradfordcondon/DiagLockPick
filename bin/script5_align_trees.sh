
mkdir ../alignmentsIndividual
mkdir ../treesIndividual

for file in *.fasta

do 

muscle3 -in $file -out ../alignmentsIndividual/$file  -maxiters 2 

done

cd alignmentsIndividual

for file in *.fasta

do 

clustalw2 -INFILE=$file -TREE 


done

for file in *.ph

do

	mv $file ../treesIndividual/

done