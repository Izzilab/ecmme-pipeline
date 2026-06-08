#!/bin/bash

WORKD="/home/petrov/Projects/ecm-evolution/pipeline/work"
WGETD="$WORKD/downloads/UniProt"
GENES="$WORKD/temp/uniprot.txt"
export WGETD

mkdir -p $WGETD

# Download UniProt features as txt
uniprot_download_txt() {
  local UniProt="${1}"
  echo "${UniProt}"
  wget \
    -t 50 \
    -c "https://www.uniprot.org/uniprot/${UniProt}.txt" \
    -O $WGETD/${UniProt}.txt
  
  # Collect the canonical isoforms  
  grep "Sequence=Displayed" $WGETD/${UniProt}.txt | awk '{print $2}' | sed 's:IsoId=::g' | sed 's:;::g' >> $WORKD/temp/UniProt.displayed.txt
}

export -f uniprot_download_txt

rm -rf $WORKD/temp/UniProt.displayed.txt
cat $GENES | while read -r a ; do
  uniprot_download_txt $a
done

rm -rf $WORKD/temp/uniprot_download_txt_empty.log
touch $WORKD/temp/uniprot_download_txt_empty.log
cd $WGETD
for i in *.txt ; do
if [ ! -s "${i}" ]; then
  echo ${i%.*} >> $WORKD/temp/uniprot_download_txt_empty.log
fi
