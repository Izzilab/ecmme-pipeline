#!/bin/bash
#set -x  # Enables debug traces
THREADS=$(nproc)

WORKD="/home/petrov/Projects/ecm-evolution/pipeline/work"
ORTHO="$WORKD/orthol"
GENES="$WORKD/output/genes-orthol.txt"

cat $GENES | while read -r g ; do
  echo ${g}.pep.fasta
  mafft-linsi --thread $THREADS $ORTHO/${g}.pep.fasta > $ORTHO/${g}.linsi.fasta
done
