#!/usr/bin/env bash
# alignAll.sh
# Initialize variable to contain the directory of the trimmed paired fastq files
fastqPath="Paired/"
# Initialize variable to contain left and right reads
leftAlign=".R1.fastq"
rightAlign=".R2.fastq"

#Initialize the output path for the aligned reads
pairedOutPath="sam/"

# Creat the output directory
mkdir -p $pairedOutPath

#Loop through left and right read fastq files in the "Paired/" directory
function alignAll {
for leftRead in $fastqPath*$leftAlign
do
    #Declare file format for left and right  read fastq files
    rightRead="${leftRead/$leftAlign/$rightAlign}"
    #Remove fastq path from the paired files and replace with nothing
    sampleName="${leftRead/$leftAlign/}"
    #Remove paired read suffix and assign to sam
    sampleName="${sampleName/$fastqPath/$pairedOutPath}"
    
    #Perform alignment function using GSNAP and output aligned reads into a .sam file
    echo nice -n19 gsnap \
        -A sam \
        -D . \
        -d AiptasiaGmapDb \
        -s AiptasiaGmapIIT.iit \
        $leftRead \
        $rightRead \
    1>$sampleName.sam
done
}
alignAll 1>alignAll.log 2>alignAll.err &


