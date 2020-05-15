#!/usr/bin/env python3

# Interleave.py

#Import Seq, SeqRecord and SeqIO
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

#Obtain fastq file for both left and right reads and parse
leftReads = SeqIO.parse("/scratch/AiptasiaMiSeq/fastq/Aip02.R1.fastq", "fastq")
rightReads = SeqIO.parse("/scratch/AiptasiaMiSeq/fastq/Aip02.R2.fastq", "fastq")

#Declare and open end file for writing
with open("interleaved.fasta","w") as output:

    #Iterate over both arrays for left and right reads
    for reads in zip(leftReads, rightReads):

        #Write reads to interleaved.fasta
        SeqIO.write(reads,output,"fastq")
