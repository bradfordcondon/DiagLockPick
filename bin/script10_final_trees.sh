
for file in *.fasta

do 

clustalw2 -INFILE=$file -TREE 


done

for file in *.ph

do

	mv $file ../treesCombined/

done