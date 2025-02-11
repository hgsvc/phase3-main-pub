#!/bin/sh

coverage=$1
fraction=$2
reads=$3
out=$4

cov=$(<$coverage)

p=$(awk -v frac=$fraction -v c=$cov 'BEGIN { printf "%.2f", frac/c }')

printf "p=%f\n" $p

[ "$fraction" -gt $cov ] && echo "Not enough coverage" && exit 1

seqtk sample $reads $p | bgzip > $out

exit 0
