#!/usr/bin/env bash
# blastPep.sh

#Align ORF amino acid sequences against the multi-organism BLAST Sequence database: SwissProt
blastp -query Trinity.fasta.transdecoder_dir/longest_orfs.pep  \
    #Declare the use of SwissProt database    
    -db swissprot  -max_target_seqs 1 \
    #Declare tabular output
    -outfmt 6 -evalue 1e-5 -num_threads 4 1> blastp.outfmt6 \
    2>blastp.err &
