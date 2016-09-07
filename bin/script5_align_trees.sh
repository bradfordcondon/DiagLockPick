
mkdir ../alignments
mkdir ../trees

for file in *.fasta

do 

muscle3 -in $file -out ../alignments/$file  -maxiters 2 

done

cd alignments

for file in *.fasta

do 

clustalw2 -INFILE=$file -TREE 


done

for file in *.ph

do

	mv $file ../trees/

done