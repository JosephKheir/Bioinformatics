#!/usr/bin/env python
#join1.py
# Define the separator
sep = "|"

# A sample kmer to use in this join
kmer = 'ATGCC'

# A sample count
count = 2

# Print the result
print(sep.join((kmer, str(count))))

