# Differential Expression Analysis of Vibrio and Menthol Data Sets versus Controls

## Methods

The goal of this script was to analyze the differential expression
levels of NGS data between Vibrio and Menthol data sets versus controls.
Three separate bioinformatics packages were used to accomplish this
which include: Salmon, tximport and DESeq2. A transcriptome index was
built using a kmer of size 25. Paired end RNA-Seq data from de-novo
assembly was quantified for relative abundance using salmon by aligning
to the transcriptome index, AipIndex. Selective alignment of the
sequence reads was accomplished using the “–validatedMappings”
parameter. In order to summarize relative abundance at the gene level, a
data frame, here called tx2gene was created in order for transcripts
from salmon quantification to be associated with the appropriate gene
IDs (Patro et al. 2017). The tximport tool functions to import
transcript abundance data from salmon and summarize to the gene level of
tx2gene (Soneson, Love, and Robinson 2015). In final, statistical
quantification for differential analysis of relative transcript
abundance data from salmon was calculated using the R package DESeq2.
DESeq2 was used to estimate the variance (i.e. equivalence or
dispersion) of NGS expression data between Vibrio and Menthol data sets
versus controls (Love, Huber, and Anders 2014). An adjusted p-value at a
significance of p \<= 0.05 was used in this study.

## Results

Based on the relative abundance of NGS transcripts, the Vibrio group
consisted of n = 58,596 samples and the Menthol group consisted of n =
15,420 samples with a total study size of 74,016 samples. Of the Vibrio
group, n = 21,588 (\~37%) of samples had a negative Log2FoldChange ratio
indicative that gene expression in this group is less relative to the
controls. Alternately, the Menthol group consisted of n = 9,252 (\~60%)
samples that had a negative Log2FoldChange ratio indicative that gene
expression in this group is less relative to controls.

## References

<div id="refs" class="references">

<div id="ref-Love">

Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15 (12): 550–50.
<https://doi.org/10.1186/s13059-014-0550-8>.

</div>

<div id="ref-Patro">

Patro, Rob, Geet Duggal, Michael I. Love, Rafael A. Irizarry, and Carl
Kingsford. 2017. “Salmon Provides Fast and Bias-Aware Quantification of
Transcript Expression.” *Nature Methods* 14 (4): 417–19.
<https://doi.org/10.1038/nmeth.4197>.

</div>

<div id="ref-Soneson">

Soneson, Charlotte, Michael I. Love, and Mark D. Robinson. 2015.
“Differential Analyses for RNA-Seq: Transcript-Level Estimates Improve
Gene-Level Inferences.” *F1000Research* 4 (December): 1521–1.
<https://www.ncbi.nlm.nih.gov/pubmed/26925227>.

</div>

</div>
