#!/usr/bin/env python
# get_one_kmer.py
seq = 'GCCGGCCCTCAGAGCAGGATGGTCCTGGATGTGTCTCCACGGACCGCAATGTGCGTTCGAAAT'
kmer_length = 6
for start in range(0, len(seq)-kmer_length):
    kmer = seq[start:start+kmer_length]
    print(kmer)
