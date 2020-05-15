#!/usr/bin/env bash
# getPathDesc.sh

#Open files for stdout and stderr
out = open('pathD.txt','w')
err = open('pathD.err','w')

wget http://rest.kegg.jp/list/pathway/ko, stdout = out, stderr = err

