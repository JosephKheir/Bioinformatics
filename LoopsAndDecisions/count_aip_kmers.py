#!/usr/bin/env python
# count_aip_kmers.py

# Import re to support regular expressions in this program
import re

# Open the fastq file for reading
read_sample = open('/scratch/AiptasiaMiSeq/fastq/Aip02.R1.fastq','r')

# Initialize a variable to contain the lines
line = ' '

# Initialize kmer length
kmer_length = 6

# Initialize the kmer dictionary
kmer_dictionary = {}

# While the line is not empty
while line:
    
    # Only read lines from the file containing sequence data
    line = read_sample.readline()
    if re.match('^[ATGCN]+$', line):

        # Get the substring (kmer) at the start of the line to length of 6 sequences
        for start in range(0, len(line) - kmer_length):
            kmer = line[start:start + kmer_length]

            # See if the kmer is in the dictionary
            if kmer in kmer_dictionary:

                # Add one to the count
                kmer_dictionary[kmer] += 1

            else:
                # It is not in the dictionary so add with a count of 1
                kmer_dictionary[kmer] = 1

# Declare the text file to write the date to
with open("aip_kmers.txt",'w') as aip_kmers:
    #  Save a tab character to output data as a tab-separated output
    t = "\t"

    # Iterate over the kmers in the dictionary
    for kmer in kmer_dictionary:

        # Get the count of how many times observed
        count = kmer_dictionary[kmer]

        # Convert the count to a string and make an array of the kmer and count
        out = [kmer, str(count)]

        #Join the elements of the out array with tabs and export to file
        aip_kmers.write(t.join(out)+"\n")
       

