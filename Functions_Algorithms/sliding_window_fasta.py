#!/usr/bin/env python3
#sliding_window_fasta.py


#Import module re and sys
import sys


#Declare a function that takes two arguments: an integer & a fasta file to
#find GC count
def sliding_window(k, string):
    """ Returns a list of k-mers from the given fasta file
    """
    kmers=[]
    end = len(string) - k + 1
    for start in range(0, end):
        kmers.append(string[start:start + k])
    #Will return kmers of size k
    return kmers


#Declare a function that takes the kmers of size k from previous function and
#count gc content
def gc_content(kmers):
    """ Returns [0,1], the fraction of GCs in the kmer
    """
    #Count the number of G's & C's
    gc = 0
    for nucleotide in kmers:
        if nucleotide in ['G', 'C']:
            gc += 1
    #Return gc concentration per kmer of length k
    return gc/len(kmers)


#Generate a code block to call the function from the CLI
if __name__=="__main__":
    #Make sure at least 2 arguments are passed
    arg_count = len(sys.argv) -1
    if arg_count < 2:
        #If less than two arguments passed, raise exception and exit function
        raise Exception("This script requires at least 2 arguments")

    
    #Declare variables to house values of first and second arguments from CLI to
    #pass to the functions
    k = int(sys.argv[1])
    kmers = {}
    #Open the fasta file passed from CLI and remove header
    with open(sys.argv[2]) as string:
        header = None
        data = ''
        for line in string:
            if line.startswith('>'):
                if header and data:
                    kmers[header] = data
                data = ''
                header = line.rstrip()
            else:
                data += line.rstrip()
        if header and data:
            kmers[header] = data
    kmers = sliding_window(k, data)
    
    
    #Calculate and print the gc concentration of kmer of length k
    for i in range(len(kmers)):
        result = gc_content(kmers[i])
        print("{}\t{:.2}".format(kmers[i], result))
