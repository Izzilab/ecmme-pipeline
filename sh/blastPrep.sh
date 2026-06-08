#!/bin/bash

WORKD="/home/petrov/Projects/ecm-evolution/pipeline/work"
export WORKD

GENES=$( cat $WORKD/output/genes-blast.txt | datamash transpose )

# Create BLAST database for each reference sequence, against which we
# will run the sequences of collected orthologues in order to determine
# best matching isoform. The reference sequence is from ncbi.
blast_db_prep() {
	local ncbiSeq=${1}
	
	makeblastdb \
		-in $WORKD/blast/hsap/$ncbiSeq/Homo_sapiens.fasta \
		-parse_seqids \
		-blastdb_version 5 \
		-title "$ncbiSeq" \
		-dbtype prot
}

# export -f the function first or use env_parallel
export -f blast_db_prep
parallel blast_db_prep ::: $GENES
