#!/bin/bash
WORKD="/home/petrov/Projects/ecm-evolution/pipeline/work"
GENES=$( cat $WORKD/output/genes-orthol.txt | datamash transpose )

SEQ="$WORKD/orthol"
GBL="$WORKD/Gblocks"
#MSA="linsi"
MSA="muscle5"

export SEQ
export GBL
export MSA

run_Gblocks(){
	local i=${1}
	Gblocks $SEQ/${i}.${MSA}.fasta -p=t -b4=5 -b5=h -s=n
	
	# Gblocks outputs to current work dir, let's move
	mv $SEQ/${i}.${MSA}.fasta-gb.txt $GBL/${i}.${MSA}.fasta-gb.txt
	
	# Prepare a "mock" sequence of only the quality. The sed insanity is from:
	# https://serverfault.com/a/391365
	grep "Gblocks          " $GBL/${i}.${MSA}.fasta-gb.txt > $GBL/$i.Gblocks.fa
	sed -i 's:Gblocks          ::g' $GBL/$i.Gblocks.fa
	sed -i ':a;N;$!ba;s/\n//g' $GBL/$i.Gblocks.fa
	sed -i '1i >Gblocks' $GBL/$i.Gblocks.fa
	echo "$i -> Gblocks"
}

mkdir -p $GBL
export -f run_Gblocks
parallel run_Gblocks ::: $GENES
