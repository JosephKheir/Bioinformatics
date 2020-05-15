## Overview

The scope of this exercise was to be able to identify variants among NGS
data compared to a reference genome. The Methods section below will
outline the purpose of each script in this module in order to accomplish
the variant calling workflow.

## Methods

getGenome.sh: The script, getGenome.sh used the function “wget” to pull
and download the GRCH38 human reference genome onto the server as a .gz
file.

getReads.sh: The next step in the variant calling workflow is to obtain
the individual NGS reads from the NA12878 reference sample. The method
invoked is the fastq-dump method that pulls the NGS data and “dumps”
each read into a separate file using the –split-files method.

trimReads.sh: The script trimReads.sh was used to quality trim paired
end reads of NGS data using the Trimmomatic method (Bolger, Lohse, and
Usadel 2014) on sequences with a minimum length of 36 bases. Parameter
“Illuminaclip” was used to trim illumina-specific sequences from the
reads. Low quality bases from both the leading and trailing end of the
reads whose scores fell below 20 were removed from the reads. A 4-base
sliding window was used to cut bases with phred quality scores that fell
below 30.

indexGenome.sh: In the indexGenome.sh script, the goal was to simply
index the reference GRCH38 genome with output in .fai file format. This
was done using the bwa index method.

alignReads.sh: The next script in the variant calling workflow is to
align NGS reads to the reference GRCH38 genome. The method invoked here
is the bwa mem algorithm which works by seeding alignments with maximal
exact matches (ie: MEM) and then extending those seeds along the
reference genome(Abuín et al. 2016). Shorter split hits were marked with
header lines passed with 8 threads as parameters during alignment.

sortSam.sh: Now that sam files from bwa have been generated and aligned
to the reference genome, they need to be sorted and create a binary sam
to bam file in order to do any meaningful evaluations. Using samtools
sort the aligned reads are sorted in genome order and saved as a .bam
file format(Li et al. 2009).

indexReads.sh: The script indexReads functions by indexing the genome
sorted BAM file. This allows for quick extraction of overlapping
alignments on particular genomic regions. The method invoked in this
script was samtools index using 8 threads.

runDeepVariant.sh: The final step in the variant calling workflow
consist of calling the observed genetic variants between NGS reads and
the reference genome. Deep variant is an analysis pipeline method using
neural network to call genetic variants (Poplin et al. 2018), and does
so in three steps in order to go from a BAM file to VCF/gVCF output
files. First, make examples: is the portion of the variant calling
script that parses through the reads and calls an alternate or “example”
allele compared to the reference. The second portion of variant calling
script is to call variants based on observed vs.  the reference genome.
The last portion of the script is the post process variants which takes
data from the examples and calls and generates genotype probabilities
and outputs as a VCF or gVCF file format (Supernat et al. 2018).

## References

<div id="refs" class="references">

<div id="ref-Abuin">

Abuín, José M., Juan C. Pichel, Tomás F. Pena, and Jorge Amigo. 2016.
“SparkBWA: Speeding up the Alignment of High-Throughput DNA Sequencing
Data.” *PloS One* 11 (5): e0155461–e0155461.
<https://doi.org/10.1371/journal.pone.0155461>.

</div>

<div id="ref-Bolger">

Bolger, Anthony M., Marc Lohse, and Bjoern Usadel. 2014. “Trimmomatic: A
Flexible Trimmer for Illumina Sequence Data.” *Bioinformatics (Oxford,
England)* 30 (15): 2114–20.
<https://doi.org/10.1093/bioinformatics/btu170>.

</div>

<div id="ref-Li">

Li, Heng, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils
Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, and 1000 Genome
Project Data Processing Subgroup. 2009. “The Sequence Alignment/Map
Format and SAMtools.” *Bioinformatics (Oxford, England)* 25 (16):
2078–9. <https://doi.org/10.1093/bioinformatics/btp352>.

</div>

<div id="ref-Poplin">

Poplin, Ryan, Pi-Chuan Chang, David Alexander, Scott Schwartz, Thomas
Colthurst, Alexander Ku, Dan Newburger, et al. 2018. “A Universal SNP
and Small-Indel Variant Caller Using Deep Neural Networks.” *Nature
Biotechnology* 36 (10): 983–87. <https://doi.org/10.1038/nbt.4235>.

</div>

<div id="ref-Supernat">

Supernat, Anna, Oskar Valdimar Vidarsson, Vidar M. Steen, and Tomasz
Stokowy. 2018. “Comparison of Three Variant Callers for Human Whole
Genome Sequencing.” *Scientific Reports* 8 (1): 17851.
<https://doi.org/10.1038/s41598-018-36177-7>.

</div>

</div>
