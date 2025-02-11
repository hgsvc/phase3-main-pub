#!/bin/bash

FULL="sort -k1,1n -k2,2n --batch-size=63 -m"
ls nonref/*.nonref.kmer.gz | split -l 40 -d - chunks
for F in chunks*
do
    BS=`cat ${F} | wc -l | cut -f 1`
    CMD="sort -k1,1n -k2,2n --batch-size=${BS} -m"
    for LINE in `cat ${F}`
    do
	CMD=${CMD}" <(zcat "${LINE}")"
    done
    CMD=${CMD}" | gzip -c > merged_${F}.gz"
    echo ${CMD} >> ${F}.sh
    chmod a+x ${F}.sh
    ./${F}.sh
    FULL=${FULL}" <(zcat merged_${F}.gz)"
done

FULL=${FULL}" | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/' | gzip -c > nonref.count.gz"
eval "${FULL}"
