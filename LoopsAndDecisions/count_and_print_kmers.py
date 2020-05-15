#!/usr/bin/env python
#count_and_print_kmers.py
seq = 'GCCGGCCCTCAGACAGGAGTGGTCCTGGATGTGGATG'
kmer_length = 6

#Initialize a k-mer dictionary
kmer_dictionary = {}

# Iterate over the positions
for start in range (0, len(seq) - kmer_length):
    #Get the substring at a specific start and end position
    kmer = seq[start:start + kmer_length]

    #see if it's in the dicitonary
    if kmer in kmer_dictionary:
        #Add one to the count
        kmer_dictionary[kmer] += 1
    else:
        #It's not in the dictionary so add with a count of 1
        kmer_dictionary[kmer] = 1

#save a tabe carachter to print output as a tab-separated
t= "\t"

#Iterate over the k-mers in the dictionary
for kmer in kmer_dictionary:
    #Get the count of how many times observed
    count = kmer_dictionary[kmer]

    #convert the count to a string and make an array of the kmer and count
    out = [kmer, str(count)]

    #Join the elements of the out array with tabs and print
    print(t.join(out))

