#!/usr/bin/env python3
#BioPython_genbank.py


#Import Entrez and SeqIO from BioPython module
from Bio import Entrez
from Bio import SeqIO


#Prove email for Entrez BioPython module
Entrez.email = "kheir.j@husky.neu.edu"


#Using BioPython Entrez, fetch from genbank sequence 515056 & J01673.1
with Entrez.efetch(
    db="nucleotide", rettype="gb", retmode="text", id="515056, J01673.1"
) as handle:
    
    
    #For the records of the above fetch sequence, print features: type, 
    #location & strand of the sequence
    for record in SeqIO.parse(handle, "gb"):
        for feature in record.features:
            print(feature.type, feature.location, feature.strand)
