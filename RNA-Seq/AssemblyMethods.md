# Transcriptome Assembly

## Joseph Kheir

## Methods

### Transcriptome Assembly

Transcriptome assembly was conducted via two methods. Method one, using Trinity (Grabherr et al. 2013) as the genome assembler, generated and assembled a transcriptome using a reference guided genome. Trinity (using the reference guided assembly method) requires a single sorted bam file for input. The package samtools (Li et al. 2009) was used to merge all sorted and indexed BAM files of RNA-Seq data into a single BAM binary alignment format file. Trinity as the transcriptome assembler operated with a maximum separation distance allowed of 1000 segments as a parameter for transcripts to span introns. Alternatively The second method of transcriptome assembly was the De-novo assembly approach. The De-novo assembly approach also used Trinity as the genome assembler, however, unlike reference guided assembly, Trinity in De-novo assembly requires comma-separated lists of files as input to assemble the transcriptome. Since short read data during initial sequencing was generated using the paired end method, comma-separated lists of quality trimmed FASTQ files (one file for all left reads, and one file for all right reads) were generated. Regardless of which method (reference guided or De-novo), Trinity FASTA output files were generated for each assembled transcriptome.

## References

Haas, B. J., Papanicolaou, A., Yassour, M., Grabherr, M., Blood, P. D., Bowden, J., … Regev, A. (2013). De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. Nature protocols, 8(8), 1494–1512. doi:10.1038/nprot.2013.084

## Citations
Grabherr, M. G., Haas, B. J., Yassour, M., Levin, J. Z., Thompson, D. A., Amit, I., … Regev, A. (2011). Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nature biotechnology, 29(7), 644–652. doi:10.1038/nbt.1883 /

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … 1000 Genome Project Data Processing Subgroup (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics (Oxford, England), 25(16), 2078–2079. doi:10.1093/bioinformatics/btp352
