#!/usr/bin/env bash
# sortAll.sh

#Initialize variable to contain directory of aligned sam files
samPath="sam/"

#Initialize .sam file extension
samExtension=".sam"

#Initialize bam  output path
bamPath="bam/"
#Declare output file type
bamExtension=".bam"
#Creat a subdirectory for the bam output path
mkdir -p $bamPath

#Loop through the aligned sam files
function sortAll {
for sortAlign in $samPath*$samExtension
    do
    #Declare files in sam subdirectory to sort
    sampleName="${sortAlign/$samPath/}"
    #Remove .sam extension and replace with nothing
    sampleName="${sampleName/$samExtension/}"
    #Sort sam files using samtools and output into bamPath as a bam file
    samtools sort $samPath$sampleName$samExtension \
        -o $bamPath$sampleName$bamExtension
    done
}
sortAll 1>sortAll.sort.log 2>sortAll.sort.err &
