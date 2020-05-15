#!/usr/bin/env python3
#find_orf.py
#import re and seq
import re
from Bio.Seq import Seq
dna=Seq("GGTCCGGGATGCCTGAATGGTACACTGGTAAGTACACTGTAAGTAAAAAAAA")
rna=dna.transcribe()
print(rna)
orf=re.search('AUG([AUGC]{3})+?(UAA|UAG|UGA)',str(rna)).group()
print(orf)

