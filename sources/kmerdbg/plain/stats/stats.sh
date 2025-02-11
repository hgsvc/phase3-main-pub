#!/bin/bash

python stats.py -k nonref1kGP.samples.ct.gz -s integrated_call_samples_v3.20130502.ALL.panel.txt > stats.1kGP.tsv

Rscript stats.R stats.1kGP.tsv
