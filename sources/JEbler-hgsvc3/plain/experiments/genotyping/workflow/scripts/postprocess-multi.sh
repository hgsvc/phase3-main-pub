#!/bin/bash
set -e

if [ $# -ne 2 ] ; then
	echo "Usage: $0 <input.vcf.gz> <output_multi_all.vcf.gz>"
	exit 1
fi

input="$1"
output_multi_all="$2"

if [ ! -f $input ] ; then
	echo "File $input not found"
	exit 1
fi

# source ~/.bashrc
# conda activate bcftools
which bcftools
bcftools view -o ${output_multi_all} -O z ${input}
tabix ${output_multi_all}
