#!/usr/bin/env bash
 
R -e "rmarkdown::render('CancerGenomics.Rmd', output_format='all')"
