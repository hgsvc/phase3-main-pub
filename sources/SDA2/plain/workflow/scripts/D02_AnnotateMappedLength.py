#!/usr/bin/env python

# Purpose: Add 'pa' tag to input bam file. The 'pa' tag indicates the percent alignment of the gene model to the gene copy in the assembly as a decimal.
# Input:
#   Arg1: bam of resolved originals (see rule D02_FilterResolvedOriginals_FiltPercentAligned in snakefile).
# Output: Input bam file with a 'pa' tag indicating alignment of gene model to assembly as decimal.
# Use: Used in rules D02 and E07 of SegDupAnnotation2

import sys
import argparse
import pathlib
import pysam

print("---- Running D02_AnnotateMappedLength.py", file=sys.stderr)
print("-- Parsing Input", file=sys.stderr)
parser = argparse.ArgumentParser(description="Add custom 'pa' tag to input bam file indicating percent alignment of gene model to gene copy in the assembly as a decimal.")
parser.add_argument("bam_filepath", type=pathlib.Path, help="Input bam file after aligning gene models to gene copies.")
args = parser.parse_args()

resOrig = pysam.AlignmentFile(args.bam_filepath)

sys.stdout.write(str(resOrig.header))

for aln in resOrig.fetch():
    astat=aln.get_cigar_stats()
    nHardClip=astat[0][5]
    nAligned = astat[0][0] + astat[0][7] + astat[0][8] # = alignment match + sequence match + sequence mismatch
    
    if (aln.query_sequence is None):
        continue
    qLen = len(aln.query_sequence) + nHardClip

    sys.stdout.write(aln.to_string() + "\tpa:f:" + str(nAligned/qLen)  + "\n")

resOrig.close()