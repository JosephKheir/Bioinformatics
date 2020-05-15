#!/usr/bin/env bash
# alignPredicted.sh

#Declare relative file path for BLAST search
blastp -query Trinity.fasta.transdecoder.pep \
    #Declare database for sequence alignment against Swissprot
    -db /blastDB/swissprot \
    #Declare output format specifications (expected value, tabular format/ length/ etc.
    -evalue 1e-10 -outfmt "6 qseqid sacc qlen slen length nident pident evalue stitle" \
    #Output data into .txt file and errors in .err file
    1>alignPredicted.txt 2>alignPredicted.err &
