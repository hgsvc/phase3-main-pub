#!/bin/bash

zcat nonref.count.gz | cut -f 1 | head -n 10000000 | gzip -c > vac.gz
Rscript vac.R vac.gz

