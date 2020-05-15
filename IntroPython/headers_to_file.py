#!/usr/bin/env python
# headers_to_file.py

# To use regular expressions, import re
import re

# Set the file path to the Drosophila genome
dmel_genome_path = '/scratch/Drosophila/dmel-all-chromosome-r6.17.fasta'

# Initialize a line counter
line_count = 0;

# Open the genome file and declare the text file to write data to
with open(dmel_genome_path) as dmel_genome, open("dmel_headers.txt",'w') as dmel_headers:

    # Creat for loop to identify genome headers
    for line in dmel_genome:

        # Check to see if the line is a header line (starts with >)
        if re.match('^>', line):

            # Increment the line count by one and limit to the first 50 lines
            line_count += 1
            if line_count < 51:

                # Write the line of output
                dmel_headers.write(line)

