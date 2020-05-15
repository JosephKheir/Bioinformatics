#!/usr/bin/env bash

# fmt_tutorial.sh


# Obtain data files
wget \
  -O "sample-metadata.tsv" \
  "https://data.qiime2.org/2020.2/tutorials/fmt/sample_metadata.tsv"

# Obtain demultiplexed sequences
wget \
  -O "fmt-tutorial-demux-1.qza" \
  "https://data.qiime2.org/2020.2/tutorials/fmt/fmt-tutorial-demux-1-10p.qza"
wget \
  -O "fmt-tutorial-demux-2.qza" \
  "https://data.qiime2.org/2020.2/tutorials/fmt/fmt-tutorial-demux-2-10p.qza"


# Perform sequence quality control using DADA2
qiime demux summarize \
  --i-data fmt-tutorial-demux-1.qza \
  --o-visualization demux-summary-1.qzv
qiime demux summarize \
  --i-data fmt-tutorial-demux-2.qza \
  --o-visualization demux-summary-2.qzv

qiime dada2 denoise-single \
  --p-trim-left 13 \
  --p-trunc-len 150 \
  --i-demultiplexed-seqs fmt-tutorial-demux-1.qza \
  --o-representative-sequences rep-seqs-1.qza \
  --o-table table-1.qza \
  --o-denoising-stats stats-1.qza
qiime dada2 denoise-single \
  --p-trim-left 13 \
  --p-trunc-len 150 \
  --i-demultiplexed-seqs fmt-tutorial-demux-2.qza \
  --o-representative-sequences rep-seqs-2.qza \
  --o-table table-2.qza \
  --o-denoising-stats stats-2.qza

# View denoising stats
qiime metadata tabulate \
  --m-input-file stats-1.qza \
  --o-visualization denoising-stats-1.qzv
qiime metadata tabulate \
  --m-input-file stats-2.qza \
  --o-visualization denoising-stats-2.qzv

# Merge Denoised stats
qiime feature-table merge \
  --i-tables table-1.qza \
  --i-tables table-2.qza \
  --o-merged-table table.qza
qiime feature-table merge-seqs \
  --i-data rep-seqs-1.qza \
  --i-data rep-seqs-2.qza \
  --o-merged-data rep-seqs.qza

# Generate merged stats summary for FeatureTable[Frequencey]
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

# Generate merged stats summary for FeatureTable[Sequences]
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

# Now analysis of metadata can be performed 
