Plink
================
Joseph Kheir
4/5/2020

# Overview

Genome wide association studies (GWAS) can bring several challenges when
it comes to computational and analytical research. The GWAS toolset
PLINK is the solution to some of these challenges. PLINK allows for
hundreds of thousands of markers that were genotyped to be analyzed and
manipulated rapidly for whole genome association studies (Purcell et al.
2007).

## Getting Started

``` bash
# Check that the input file is intact
plink --file hapmap1
```

## Making a binary PED file

``` bash
# Make a binary PED file
plink --file hapmap1 --make-bed --out hapmap1

# Making a binary PED file with selecting high genotyping
plink --file hapmap1 --make-bed --mind 0.05 --out highgeno
```

``` r
# When using the --make-bed option, the threshold filters for missing rates and allele frequency were automatically set to exclude nobody
library(knitr)
table <- read.table("hapmap1.fam", sep = "\t")
knitr::kable(head(table), "markdown")
```

| V1               |
| :--------------- |
| HCB181 1 0 0 1 1 |
| HCB182 1 0 0 1 1 |
| HCB183 1 0 0 1 2 |
| HCB184 1 0 0 1 1 |
| HCB185 1 0 0 1 1 |
| HCB186 1 0 0 1 1 |

``` r
table1 <- read.table("hapmap1.bim", sep = "\t")
knitr::kable(head(table1), "markdown")
```

| V1 | V2        | V3 | V4 | V5 | V6 |
| -: | :-------- | -: | -: | -: | -: |
|  1 | rs6681049 |  0 |  1 |  1 |  2 |
|  1 | rs4074137 |  0 |  2 |  1 |  2 |
|  1 | rs7540009 |  0 |  3 |  0 |  2 |
|  1 | rs1891905 |  0 |  4 |  1 |  2 |
|  1 | rs9729550 |  0 |  5 |  1 |  2 |
|  1 | rs3813196 |  0 |  6 |  1 |  2 |

## Summary statistics: missing rates

``` bash
# Generate summary statistics on rates of missing data
plink --bfile hapmap1 --missing --out miss_stat

# Use "-more" command to view data
more miss_stat.lmiss
more miss_stat.imiss

# Analyze the data by chromosome
plink --bfile hapmap1 --chr 1 --out res1 --missing
plink --bfile hapmap1 --chr 2 --out res2 --missing
```

``` r
# summary statistics on rates of missing data in the file, using the --missing option
library(knitr)
table <- read.table("miss_stat.lmiss", sep = "\t")
knitr::kable(head(table), "markdown")
```

| V1                              |
| :------------------------------ |
| CHR SNP N\_MISS N\_GENO F\_MISS |
| 1 rs6681049 0 89 0              |
| 1 rs4074137 0 89 0              |
| 1 rs7540009 0 89 0              |
| 1 rs1891905 0 89 0              |
| 1 rs9729550 0 89 0              |

``` r
table1 <- read.table("miss_stat.imiss", sep = "\t")
knitr::kable(head(table1), "markdown")
```

| V1                                          |
| :------------------------------------------ |
| FID IID MISS\_PHENO N\_MISS N\_GENO F\_MISS |
| HCB181 1 N 671 83534 0.008033               |
| HCB182 1 N 1156 83534 0.01384               |
| HCB183 1 N 498 83534 0.005962               |
| HCB184 1 N 412 83534 0.004932               |
| HCB185 1 N 329 83534 0.003939               |

## Summary statistics: allele frequencies

``` bash
# Performing summary statistics on alelle frequencies
plink --bfile hapmap1 --freq --out freq_stat

# Use "within" to perform stratified analysis
plink --bfile hapmap1 --freq --within pop.phe --out freq_stat
more freq_stat.frq.strat

# Use -snp to get statistics for specific SNP
plink --bfile hapmap1 --snp rs1891905 --freq --within pop.phe --out snp1_frq_stat
```

``` r
# frequency analysis (and the missingness analysis) stratified by a categorical, cluster variable
library(knitr)
table <- read.table("freq_stat.frq.strat", sep = "\t")
knitr::kable(head(table), "markdown")
```

| V1                                 |
| :--------------------------------- |
| CHR SNP CLST A1 A2 MAF MAC NCHROBS |
| 1 rs6681049 1 1 2 0.2333 21 90     |
| 1 rs6681049 2 1 2 0.1932 17 88     |
| 1 rs4074137 1 1 2 0.1 9 90         |
| 1 rs4074137 2 1 2 0.05682 5 88     |
| 1 rs7540009 1 0 2 0 0 90           |

## Basic association analysis

``` bash
# Perform basic association analysis
plink --bfile hapmap1 --assoc --out as1

# Sort association list and print first 10
sort --key=7 -nr as1.assoc | head

# A sorted association list with a range of segnificant values
plink --bfile hapmap1 --assoc --adjust --out as2
more as2.assoc.adjusted

# Test for frequencies between groups
plink --bfile hapmap1 --pheno pop.phe --assoc --adjust --out as3
```

``` r
# Here we see that the simulated disease variant rs2222162 is actually the second most significant SNP in the list, with a large difference in allele frequencies of 0.28 in cases versus 0.62 in controls
library(knitr)
table <- read.table("as2.assoc.adjusted", sep = "\t")
knitr::kable(head(table), "markdown")
```

| V1                                                                     |
| :--------------------------------------------------------------------- |
| CHR SNP UNADJ GC BONF HOLM SIDAK\_SS SIDAK\_SD FDR\_BH FDR\_BY         |
| 13 rs9585021 5.586e-06 4.994e-05 0.3839 0.3839 0.3188 0.3188 0.09719 1 |
| 2 rs2222162 5.918e-06 5.232e-05 0.4068 0.4067 0.3342 0.3342 0.09719 1  |
| 9 rs10810856 7.723e-06 6.483e-05 0.5308 0.5308 0.4118 0.4118 0.09719 1 |
| 2 rs4675607 8.05e-06 6.703e-05 0.5533 0.5533 0.4249 0.4249 0.09719 1   |
| 2 rs4673349 8.485e-06 6.994e-05 0.5832 0.5831 0.4419 0.4419 0.09719 1  |

## Genotypic and other association models

``` bash
# Association statistics using 2-by-3 genotype table & standard allelic test using mode1
plink --bfile hapmap1 --model --snp rs2222162 --out mod1
# Declare minimum number per cell
plink --bfile hapmap1 --model --cell 0 --snp rs2222162 --out mod2
```

## Stratification analysis

``` bash
# Cluster analysis on basis of genetic identity
plink --bfile hapmap1 --cluster --mc 2 --ppc 0.05 --out str1
more str1.cluster1
```

``` r
library(knitr)
table <- read.table("str1.cluster1", sep = "\t")
knitr::kable(head(table), "markdown")
```

| V1    | V2                  |
| :---- | :------------------ |
| SOL-0 | HCB181\_1 JPT260\_1 |
| SOL-1 | HCB182\_1 HCB225\_1 |
| SOL-2 | HCB183\_1 HCB194\_1 |
| SOL-3 | HCB184\_1 HCB202\_1 |
| SOL-4 | HCB185\_1 HCB217\_1 |
| SOL-5 | HCB186\_1 HCB201\_1 |

## Association analysis accounting for clusters

``` bash
plink --bfile hapmap1 --mh --within str1.cluster2 --adjust --out aac1
more aac1.cmh.adjusted
# Pair up the most similar individuals
plink --bfile hapmap1 --cluster --cc --ppc 0.01 --out version2
# Alternate cluster scheme
plink --bfile hapmap1 --mh --within version2.cluster2 --adjust --out aac2
# Specify the number of clusters wanted
plink --bfile hapmap1 --cluster --K 2 --out version3
# External clustering
plink --bfile hapmap1 --mh --within pop.phe --adjust --out aac3

# Visualize substrate by making patrix of pairwise IBS distances
plink --bfile hapmap1 --cluster --matrix --out ibd_view
```

## Then in R, perform the following analysis on the above association analysis

``` r
m <- as.matrix(read.table("ibd_view.mibs"))
mds <- cmdscale(as.dist(1-m))
k <- c( rep("green",45) , rep("blue",44) )
plot(mds,pch=20,col=k)
```

![](runPlink_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

## Quantitative trait association analysis

``` bash
# By method of least squares
plink --bfile hapmap1 --assoc --pheno qt.phe --out quant1
```

``` r
library(knitr)
table <- read.table("quant1.qassoc", sep = "\t")
knitr::kable(head(table), "markdown") 
```

| V1                                                      |
| :------------------------------------------------------ |
| CHR SNP BP NMISS BETA SE R2 T P                         |
| 1 rs6681049 1 89 -0.2266 0.3626 0.004469 -0.6249 0.5336 |
| 1 rs4074137 2 89 -0.2949 0.6005 0.002765 -0.4911 0.6246 |
| 1 rs7540009 3 89 NA NA NA NA NA                         |
| 1 rs1891905 4 89 -0.1053 0.3165 0.001272 -0.3328 0.7401 |
| 1 rs9729550 5 89 0.5402 0.4616 0.0155 1.17 0.2451       |

``` bash
# By method of purmutation
plink --bfile hapmap1 --assoc --pheno qt.phe --perm --within str1.cluster2 --out quant2
```

``` bash
# Allow for multiple testing
plink --bfile hapmap1 --assoc --pheno qt.phe --mperm 1000 --within str1.cluster2 --out quant3
plink --bfile hapmap1 --pheno qt.phe --gxe --covar pop.phe --snp rs2222162 --out quant3
```

## Extracting a SNP of interest

``` bash
plink --bfile hapmap1 --snp rs2222162 --recodeAD --out rec_snp1
```

``` r
d <- read.table("rec_snp1.raw" , header=T)
summary(glm(PHENOTYPE-1 ~ rs2222162_1, data=d, family="binomial"))
```

## Citations

<div id="refs" class="references">

<div id="ref-Purcell">

Purcell, Shaun, Benjamin Neale, Kathe Todd-Brown, Lori Thomas, Manuel A.
R. Ferreira, David Bender, Julian Maller, et al. 2007. “PLINK: A Tool
Set for Whole-Genome Association and Population-Based Linkage Analyses.”
*American Journal of Human Genetics* 81 (3): 559–75.
<https://doi.org/10.1086/519795>.

</div>

</div>
