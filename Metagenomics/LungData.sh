#!/usr/bin/env bash
R -e "rmarkdown::render('lungData.Rmd', output_format='all')"
