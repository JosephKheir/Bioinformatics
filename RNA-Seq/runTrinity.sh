#!/usr/bin/env bash
# runTrinity.sh

#Declare path of trinity assembler
nice -n19 /usr/local/programs/trinityrnaseq-Trinity-v2.8.4/Trinity \

#Declare for Trinity to use the genome guided assembly method using the merged bam file
--genome_guided_bam bam/AipAll.bam \
#Declare maximum separation distance that Trinity will allow for segments of transcripts to span intron sites
--genome_guided_max_intron 10000 \
#Declare maximum memory Trinity is allowed to use for the assembly process
--max_memory 10G --CPU 4 \
1>trinity.log 2>trinity.err &
