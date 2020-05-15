#!/usr/bin/env python3
#BioPython_seqio.py


#Import SeqIO, Seq, SeqRecord and sys from respective modules
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys


#Declare variable to house the opened fasta file (passed in the CLI
#containing the sequence
sequence = SeqIO.parse(sys.argv[1], "fasta")


#Declare and open the file name passed in the CLI to house the reverse
#complement of the sequence from above
with open(sys.argv[2], "w"):
    for seq in sequence:
        reverse = seq.reverse_complement()
        #Write the reverse complement to the new file
        SeqIO.write(reverse, sys.argv[2], "fasta")


#Check to make sure at least two arguments are given at the CLI
if __name__=="__main__":
    arg_count = len(sys.argv) -1
    if arg_count < 2:
        raise Exception("This script requires at least 2 arguments")

