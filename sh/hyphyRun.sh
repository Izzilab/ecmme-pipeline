#!/bin/bash

WORKD="/home/petrov/Projects/ecm-evolution/pipeline/work"
ORTHO="$WORKD/orthol"
TREES="$WORKD/trees"
GENES="$WORKD/output/genes-orthol.txt"

HYPHY="/var/tmp/hyphy"

HPBIN="HYPHYMPI"
#HPBIN="hyphy"

mkdir -p $HYPHY
cd $HYPHY

cat $GENES | while read -r gene ; do
	mkdir -p $HYPHY/${gene}/{meme,fubar}
	cd $HYPHY/${gene}/meme
	$HPBIN meme --alignment $ORTHO/${gene}.codon.fasta --tree $TREES/${gene}.tre --output ${gene}.meme.json
	
	cd $HYPHY/${gene}/fubar
	$HPBIN fubar --alignment $ORTHO/${gene}.codon.fasta --tree $TREES/${gene}.tre --output ${gene}.fubar.json
	
	# log when done
	echo "$(date) $gene" >> $WORKD/temp/hyphy.log 
	cd $HYPHY
done
