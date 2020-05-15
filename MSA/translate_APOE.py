#!/usr/bin/env python3

# translate_APOE.py


# Import appropriate BioPython and re modules
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import re


# Declare a function that translates the given sequences
def translate():
    """ A function that accepts a fasta file input of multiple
    sequences and returns a new fasta file with the translated
    protein. 
    Args: APOE_refseq_transcript.fasta: fasta file of multiple 
    sequences
    """

    with open("APOE_refseq_transcript.fasta") as myfile:
        with open("apoe_aa.fasta", "w") as output:
            # Parese through APOE sequences
            for record in SeqIO.parse(myfile, "fasta"):
                rna = record.seq
                # Translate the APOE sequences to protein sequences
                protein = str(rna.translate())
                out = (SeqRecord(Seq(protein, IUPAC.protein), 
                    id=record.id, description="", name=""))
                
                # Write the protein sequences to output file
                SeqIO.write(out, output, "fasta")


# Call and run the function
translate()
