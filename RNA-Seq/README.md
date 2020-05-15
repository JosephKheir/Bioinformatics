# Short Read Alignment

## Joseph Kheir

## Methods

### Quality control and read alignment
QC sequencing runs were used in the analysis of a total of 24 anemones. Quality trimming of sequences was performed using Trimmomatic with paired-end reads on sequences with a minimum length of 36 bases. Parameter "ILLUMINACLIP" was used to trim Illumina-specific sequences from the reads. Additionally, low quality bases whose phred quality scores fell below 20 were removed from both the leading and trailing end of the reads. A 4-base sliding window was used to cut bases with phred quality scores below 30. A GMAP database and indexed intron splice sites were built and used by GSNAP to align trimmed reads. The SAM file format was used for each paired aligned read output. Alignments were sorted and indexed using samtools with BAM version as the output file parameter.

## References
http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

http://research-pub.gene.com/gmap/src/README

https://samtools.github.io/hts-specs/SAMv1.pdf 

## Citations
https://onesearch.library.northeastern.edu/primo-explore/fulldisplay?docid=TN_oxford10.1093%2Fbioinformatics%2Fbtu170&context=PC&vid=NU&lang=en_US&search_scope=default_scope&adaptor=primo_central_multiple_fe&tab=default_tab&query=any,contains,trimmomatic&offset=0
