#!/usr/bin/env python3
#startseq.py
#import re
import re
#initialize RNA seq
rna='GGUCCGGGAUGCCUGAAUGGUACACUGGUAAGUACACUGUAAGUAAAAAA'
#search for aug followed by character

orf=re.search('AUG[AUGC]{3})+?(UAA|UAG|UGA)',str(rna)).group()

#print the matching string
print(orf)
