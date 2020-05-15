#!/usr/bin/env python
#join.py
#Define the separator character
sep = "|"

# Asmpla kmer to use in this join
kmer= 'ATGCC'

#Asample count
count = 2

#convert the count to a string
str_count = str(count)

#Put the kmer and count in an array
kmer_count = [kmer, str_count]

# Join the elements of the array with the separator and save in a string
joined_string = sep.join(kmer_count)

#print the result
print(joined_string)

