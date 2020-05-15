#!/usr/bin/env python3

#find_dmel_orf.py

#Import re and Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re

#read fasta file
for record in SeqIO.parse("/scratch/Drosophila/dmel-all-chromosome-r6.17.fasta", "fasta"):
        if re.match("^\d{1}\D*$", record.id):
            dna=record.seq
            rna=dna.transcribe()
            seq=Seq(str(rna))
            orf=re.search('AUG([AUGC]{3})+?(UAA|UAG|UGA)',str(seq)).group()
            protein=Seq(orf).translate()
            print(protein)
