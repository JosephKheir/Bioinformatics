#!/usr/bin/env python
# count_kmers.py
seq = 'GCCGGCCCTCAGACAGGAGTGGTCCTGGATGTGGATG'
kmer_length = 6

# Initialize a k-mer dictionar
kmer_dictionary = {}
# Iterate over the positions
for start in range(0, len(seq) - kmer_length):
    # Get the substring at a specific start and end position
    kmer = seq[start:start+ kmer_length]

    #see if it's in the dicitonary
    if kmer in kmer_dictionary:
        #Add one to the count
        kmer_dictionary[kmer] += 1
    else:
        #It's not in the dictionary so add with a count of 1
        kmer_dictionary[kmer] = 1

#Print the number of keys in the dictionary
print(len(kmer_dictionary))

