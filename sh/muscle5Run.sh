#!/bin/bash

WORKD="/home/petrov/Projects/ecm-evolution/pipeline/work"
ORTHO="$WORKD/orthol"
GENES="$WORKD/output/genes-orthol.txt"

cat $GENES | while read -r g ; do
  echo ${g}.pep.fasta
  muscle5 -align $ORTHO/${g}.pep.fasta -output $ORTHO/${g}.muscle5.fasta
done
