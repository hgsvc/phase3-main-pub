#!/bin/sh

set -e

bamfile=$1
outfile=$2

length=$(python3 workflow/scripts/compute-genome-length.py $bamfile)
printf "computed genome length: %d\n" $length
samtools depth -a $bamfile | awk -v len=$length '{sum+=$3} END {print sum/len}' > $outfile
