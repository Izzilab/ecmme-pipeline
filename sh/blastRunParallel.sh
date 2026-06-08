#!/bin/bash

# Run reciprocal BLASTP against reference protein and sort by bit score.
# Then sort by "bit score" (12th column, --key 12). Output format info:
# http://www.metagenomics.wiki/tools/blast/blastn-output-format-6

WORKD="/home/petrov/Projects/ecm-evolution/pipeline/work"
export WORKD

GENES="$WORKD/output/genes-blast.txt"

# This is the old order of BLASTP results:     6 qseqid sacc pident ppos length mismatch gapopen gaps qstart qend sstart send evalue bitscore"
# This is the default order of BLASTP results: qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
# This is current (6 = Tabular; 10 = Comma-separated values)
#    	    qseqid means Query Seq-id
#    	       qgi means Query GI
#    	      qacc means Query accession
#    	   qaccver means Query accession.version
#    	      qlen means Query sequence length
#    	    sseqid means Subject Seq-id
#    	 sallseqid means All subject Seq-id(s), separated by a ';'
#    	       sgi means Subject GI
#    	    sallgi means All subject GIs
#    	      sacc means Subject accession
#    	   saccver means Subject accession.version
#    	   sallacc means All subject accessions
#    	      slen means Subject sequence length
#    	    qstart means Start of alignment in query
#    	      qend means End of alignment in query
#    	    sstart means Start of alignment in subject
#    	      send means End of alignment in subject
#    	      qseq means Aligned part of query sequence
#    	      sseq means Aligned part of subject sequence
#    	    evalue means Expect value
#    	  bitscore means Bit score
#    	     score means Raw score
#    	    length means Alignment length
#    	    pident means Percentage of identical matches
#    	    nident means Number of identical matches
#    	  mismatch means Number of mismatches
#    	  positive means Number of positive-scoring matches
#    	   gapopen means Number of gap openings
#    	      gaps means Total number of gaps
#    	      ppos means Percentage of positive-scoring matches
#    	    frames means Query and subject frames separated by a '/'
#    	    qframe means Query frame
#    	    sframe means Subject frame
#    	      btop means Blast traceback operations (BTOP)
#    	    staxid means Subject Taxonomy ID
#    	  ssciname means Subject Scientific Name
#    	  scomname means Subject Common Name
#    	sblastname means Subject Blast Name
#    	 sskingdom means Subject Super Kingdom
#    	   staxids means unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)
#    	 sscinames means unique Subject Scientific Name(s), separated by a ';'
#    	 scomnames means unique Subject Common Name(s), separated by a ';'
#    	sblastnames means unique Subject Blast Name(s), separated by a ';' (in alphabetical order)
#    	sskingdoms means unique Subject Super Kingdom(s), separated by a ';' (in alphabetical order) 
#    	    stitle means Subject Title
#    	salltitles means All Subject Title(s), separated by a '<>'
#    	   sstrand means Subject Strand
#    	     qcovs means Query Coverage Per Subject
#    	   qcovhsp means Query Coverage Per HSP
#    	    qcovus means Query Coverage Per Unique Subject (blastn only)
#	 

# 		-gapopen 10 -gapextend 1 \
   
blastRun(){
	local i=${1}
	local j=${2}
	blastp \
		-query $j \
		-db $WORKD/blast/hsap/$i/Homo_sapiens.fasta \
		-out $WORKD/blast/hsap/$i/results/${j%.fasta}.csv \
		-outfmt="10 qseqid sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore score ppos qlen slen gaps qcovs"
	echo "$i -> $j"
}
export -f blastRun

# These will be screened
cd $WORKD/blast/nohsap

cat $GENES | while read -r i ; do
	mkdir -p $WORKD/blast/hsap/$i/results
	cd $WORKD/blast/nohsap/$i
	parallel "blastRun $i" ::: *
	#cd ..
done
