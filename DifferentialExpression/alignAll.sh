#!/usr/bin/env bash
# alignAll1.sh


# Declare the output directory name
outDir='quant/'


# Declare the file path where samples are stored
fastqPath="/scratch/SampleDataFiles/Paired/"


# Declare sample prefix with wild-card
prefix="Aip*"


# Declare variable for left and right reads
leftRead=".R1.paired.fastq"
rightRead=".R2.paired.fastq"


# Declare a function to calculate salmon statistics
function align {
for file in $fastqPath*$leftRead
do
    # Remove trailing suffix from sample name
    suffix="${file/$fastqPath/}"
    sampleName="${suffix/$leftRead/}"
    # Run salmon aligning against AiIndex and output 
    # results to output directory
    salmon quant -l IU \
        -1 $fastqPath$sampleName.R1.paired.fastq \
        -2 $fastqPath$sampleName.R2.paired.fastq \
        -i AipIndex \
        --validateMappings \
        -o $outDir$sampleName
done
}

align 1>align.log 2>align.err &
