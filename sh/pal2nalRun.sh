#!/bin/bash
WORKD="/home/petrov/Projects/ecm-evolution/pipeline/work"
GENES=$( cat $WORKD/output/genes-orthol.txt | datamash transpose )
#MSA="linsi"
MSA="muscle5"
export MSA

run_pal2nal(){
	local i=${1}
	pal2nal.pl $i.$MSA.fasta $i.cds.fasta -output fasta > $i.codon.fasta
	echo "$i -> codon"
}

cd $WORKD/orthol
export -f run_pal2nal
parallel run_pal2nal ::: $GENES
