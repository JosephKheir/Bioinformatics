#!/usr/bin/env python3
#sliding_window_fasta.py


#Import module SeqIO/ SeqRecord from BioPython & sys
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


#Declare variables to house input from CLI
k = int(sys.argv[1])
fasta_file = sys.argv[2]


#Declare a function to return kmers of size k from fasta file passed
#in the CLI
def sliding_window(k, fasta_file):
    """ A function that returns a kmer of size k from the provided 
    fasta file
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        kmers = []
        end = len(record.seq) - k + 1
        for start in range(0, end):
            kmers.append(record.seq[start:start + k])
        return(kmers)


#Declare a function to count GC content from kmers of size k from 
#previous function
def gc_content(kmers):
    """ Returns the GC content of a kmer of size k
    """
    gc = 0
    for nucleotide in kmers:
        if nucleotide in ['G', 'C']:
            gc += 1
    return gc/len(kmers)


#Check to make sure at least two arguments are passed in the CLI
if __name__=="__main__":
    arg_count = len(sys.argv) - 1
    if arg_count < 2:
        raise Exception("This script requires at least two arguments")
    
    
    #Print a tab separated list of the kmers of size k and its corresponding 
    #GC content
    kmers = {}
    kmers = sliding_window(k, fasta_file)
    for i in range(len(kmers)):
        gc = gc_content(kmers[i])
        print("{}\t{:.2}".format(kmers[i], gc))
