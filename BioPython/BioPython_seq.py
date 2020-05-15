#!/usr/bin/env python3
#BioPython_seq.py


#Import BioPython modules
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO


#Generate variable to hold sequence record and change to uppercase
seq = Seq("aaaatgggggggggggccccgtt", generic_dna)
seq = seq.upper()


#Declare a record as a SeqRecord for the above sequence and provide
#id, and description
record = SeqRecord(seq, id = "#1234", description = "example 1")


#Write the record of the sequence to a genbank file
SeqIO.write(record, "BioPython_seq.gb", "genbank")
