#!/bin/bash

WORKD="/home/petrov/Projects/ecm-evolution/pipeline/work"
export WORKD

# Genes
GENES=$( cat $WORKD/output/genes.txt | datamash transpose )

# dump everything here
mkdir -p "$WORKD/genes/tmp"

# A function to download summary info
download_info(){
  local gene="${1}"
  
  datasets summary gene gene-id $gene --ortholog 'placentalia' --report product --as-json-lines > $WORKD/genes/tmp/$gene.jsonl
  echo "${gene}.jsonl"

  ### TODO: learn the fields!!!
  # https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/dataformat/tsv/dataformat_tsv_gene-product/
  dataformat tsv gene-product --inputfile $WORKD/genes/tmp/$gene.jsonl --fields gene-id,tax-id,tax-name,symbol,transcript-accession,transcript-length,transcript-protein-accession,transcript-protein-length,transcript-protein-isoform > $WORKD/genes/${gene}.tsv
  echo "${gene}.tsv"
}

# A function to download sequences
download_seq() {
  local gene="${1}"

  # sequences
  # https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/gene/
  # no need to download protein?
  datasets download gene gene-id $gene --ortholog placentalia --include cds --filename "$WORKD/genes/tmp/$gene.zip"
  #datasets download gene gene-id $gene --ortholog placentalia --include cds,protein --filename "$WORKD/genes/tmp/$gene.zip"
  echo "${gene}.zip"
  
  mkdir -p $WORKD/genes/tmp/$gene
  bsdtar -xf $WORKD/genes/tmp/$gene.zip -C $WORKD/genes/tmp/$gene
  cp $WORKD/genes/tmp/$gene/ncbi_dataset/data/cds.fna $WORKD/genes/$gene.fna
  #cp $WORKD/genes/tmp/$gene/ncbi_dataset/data/protein.faa $WORKD/genes/$gene.faa
}

# export -f the function first or use env_parallel
export -f download_info
export -f download_seq

parallel download_info ::: $GENES
parallel download_seq ::: $GENES

# A function to check for empty summary files
log_empty_json(){
  rm -rf $WORKD/temp/jsonl_empty.log
  touch $WORKD/temp/jsonl_empty.log
  
  cd $WORKD/genes/tmp
  for i in *.jsonl ; do
    if [ ! -s "${i}" ]; then
      echo ${i%.*} >> $WORKD/temp/jsonl_empty.log
    fi
  done
}

# A function to check for empty dirs
log_empty_dirs(){
  rm -rf $WORKD/temp/dirs_empty.log
  touch $WORKD/temp/dirs_empty.log
  
  cd $WORKD/genes/tmp
  find ./ -maxdepth 1 -type d -empty -exec echo {} >> $WORKD/temp/dirs_empty.log \;
  sed -i 's:\.\/::g' $WORKD/temp/dirs_empty.log
}

export -f log_empty_json
export -f log_empty_dirs

log_empty_json
log_empty_dirs

# Run until everything is downloaded and the log files are thus empty
GENES=$( cat $WORKD/temp/jsonl_empty.log | datamash transpose )
until [ ! -s $WORKD/temp/jsonl_empty.log ]; do
  parallel download_info ::: $GENES
  log_empty_json
done

GENES=$( cat $WORKD/temp/dirs_empty.log | datamash transpose )
until [ ! -s $WORKD/temp/dirs_empty.log ]; do
  parallel download_seq ::: $GENES
  log_empty_dirs
done
