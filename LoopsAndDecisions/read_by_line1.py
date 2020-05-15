#!/usr/bin/env python
#read_by_line1.py
#Open a sample fastq file for reading
#Import re to support regular expressions in this program. 
import re
read_sample = open('/scratch/SampleDataFiles/Sample.R1.fastq','r')

#Initialize a variable to containe the lines
line = ' '
#While line is not empty
while line:
    #Read one line from the file
    line = read_sample.readline()
    if re.match('^[ATGCN]+$',line):
        #print the line
        print(line)

#close the file
read_sample.close()

