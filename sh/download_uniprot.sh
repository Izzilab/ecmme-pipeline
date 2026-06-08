#!/bin/bash

WORKD="/home/petrov/Projects/ecm-evolution/pipeline/work"
WGETD="$WORKD/downloads/UniProt"
GENES="$WORKD/temp/UniProt.canonical.txt"
export WGETD

mkdir -p $WGETD

uniprot_download() {
  local UniProt="${1}"
  echo "${UniProt}"
  wget \
    -t 50 \
    -c "https://www.uniprot.org/uniprot/${UniProt}.fasta" \
    -O $WGETD/${UniProt}.fa
}

export -f uniprot_download

cat $GENES | while read -r a ; do
  uniprot_download $a
done

rm -rf $WORKD/temp/uniprot_download_fasta_empty.log
touch $WORKD/temp/uniprot_download_fasta_empty.log
cd $WGETD
for i in *.fa ; do
if [ ! -s "${i}" ]; then
  echo ${i%.*} >> $WORKD/temp/uniprot_download_fasta_empty.log
fi
