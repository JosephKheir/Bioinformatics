#!/usr/bin/env bash
# indexAll.sh
#Initialize the bam subdirectory for files to index
bamPath="bam/"
bamExtension=".bam"

#Loop through sorted bam files
function indexAll {
    for indexBam in $bamPath*$bamExtension
    do
        #Index sorted bam files using samtools
        samtools index $indexBam
    done
}
indexAll 1>indexAll.log 2>indexAll.err &
