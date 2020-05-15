#!/usr/bin/env python3
# get_path.py

import subprocess

#Open the KO mappings
ko = open('ko.txt')
ko_ids = {}

print(ko_ids)

#Open files for output and errors
out = open('path.txt','w')
err = open('path.err','w')

#Iterate over lines of KO ID's
for ko_line in ko:
    
    #Remove line terminator
    ko_line = ko_line.rstrip()

    print(ko_line)

    #Split line on whitespace
    fields = ko_line.split()
    
    if len(fields) > 1:
        ko = fields[1]
        ko_ids[ko] = 1
    
for ko_id in ko_ids:
    print(ko_id)
    result = subprocess.call("curl http://rest.kegg.jp/link/pathway/ko/" + ko_id, stdout=out, stderr=err, shell=True)
