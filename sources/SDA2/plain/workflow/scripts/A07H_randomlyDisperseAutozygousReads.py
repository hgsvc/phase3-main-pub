#!/usr/bin/env python

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 02/20/24

# Purpose: Randomly split autozygous reads. 
# Assumption: input bams are in lexicographic order, and PG tag already exists in input bam headers
# Input: 
# Output: 

import sys
import argparse
import pathlib
import pysam
from random import randint

# Parse Input
print("---- Running A07H_randomlyDisperseAutozygousReads.py", file=sys.stderr)
print("-- Parsing Input", file=sys.stderr)
parser = argparse.ArgumentParser(description="Randomly split aligned autozygous reads between hap1 and hap2 phased reads.")
parser.add_argument("hap1_only_bam_in", type=pathlib.Path, help="Path to bam of reads phased to hap1. Data will be copied to hap1_bam_out.")
parser.add_argument("hap2_only_bam_in", type=pathlib.Path, help="Path to bam of reads phased to hap2. Data will be copied to hap2_bam_out.")
parser.add_argument("hap1_auto_bam_in", type=pathlib.Path, help="Path to bam of autozygous reads aligned to hap1.")
parser.add_argument("hap2_auto_bam_in", type=pathlib.Path, help="Path to bam of autozygous reads aligned to hap2. Must have same exact reads as hap1_auto_bam_in.")
parser.add_argument("hap1_bam_out", type=pathlib.Path, help="Output bam filepath of reads phased to hap1 plus statistically half of purely autozygous input reads. (Output is unsorted.)")
parser.add_argument("hap2_bam_out", type=pathlib.Path, help="Output bam filepath of reads phased to hap2 plus statistically (other) half of purely autozygous input reads. (Output is unsorted.)")
args = parser.parse_args()

### Open Files and Generate Headers
print("-- Open Files and Generate Headers", file=sys.stderr)
## Open Input Files
hap1_bam_in = pysam.AlignmentFile(args.hap1_only_bam_in, "rb")
hap2_bam_in = pysam.AlignmentFile(args.hap2_only_bam_in, "rb")
hap1_auto_bam_in = pysam.AlignmentFile(args.hap1_auto_bam_in, "rb")
hap2_auto_bam_in = pysam.AlignmentFile(args.hap2_auto_bam_in, "rb")
## Create Header Dictionaries
hap1_header_dict_updated = hap1_bam_in.header.as_dict()
hap2_header_dict_updated = hap2_bam_in.header.as_dict()
hap1_header_dict_updated["HD"]["SO"]="unsorted"
hap2_header_dict_updated["HD"]["SO"]="unsorted"
hap1_bam_last_PG_ID = hap1_header_dict_updated["PG"][-1]["ID"]
hap2_bam_last_PG_ID = hap2_header_dict_updated["PG"][-1]["ID"]
hap1_new_PG_entry = {'ID':'A07H_randomlyDisperseAutozygousReads', 'PN':'A07H_randomlyDisperseAutozygousReads.py', 'PP':hap1_bam_last_PG_ID,'CL':" ".join(sys.argv)}
hap2_new_PG_entry = {'ID':'A07H_randomlyDisperseAutozygousReads', 'PN':'A07H_randomlyDisperseAutozygousReads.py', 'PP':hap2_bam_last_PG_ID,'CL':" ".join(sys.argv)}
hap1_header_dict_updated["PG"].append(hap1_new_PG_entry)
hap2_header_dict_updated["PG"].append(hap2_new_PG_entry)
## Open Output Files
hap1_bam_out = pysam.AlignmentFile(args.hap1_bam_out, 'wb', header=hap1_header_dict_updated)
hap2_bam_out = pysam.AlignmentFile(args.hap2_bam_out, 'wb', header=hap2_header_dict_updated)

## Iterate through bams
print("-- Iterate through bams", file=sys.stderr)
hap1_iter = hap1_bam_in.fetch(until_eof=True)
hap2_iter = hap2_bam_in.fetch(until_eof=True)
hap1_auto_iter = hap1_auto_bam_in.fetch(until_eof=True)
hap2_auto_iter = hap2_auto_bam_in.fetch(until_eof=True)

hap1_read = next(hap1_iter, None)
hap2_read = next(hap2_iter, None)
hap1_auto_read = next(hap1_auto_iter, None)
hap2_auto_read = next(hap2_auto_iter, None)

while hap1_read:
    hap1_bam_out.write(hap1_read)
    hap1_read = next(hap1_iter, None)

while hap2_read:
    hap2_bam_out.write(hap2_read)
    hap2_read = next(hap2_iter, None)

while hap1_auto_read and hap2_auto_read:
    if randint(0,1):
        hap1_bam_out.write(hap1_auto_read)
    else:
        hap2_bam_out.write(hap2_auto_read)
    hap1_auto_read = next(hap1_auto_iter, None)
    hap2_auto_read = next(hap2_auto_iter, None)

### Close Files
print("-- Close Files", file=sys.stderr)
hap1_bam_in.close()
hap2_bam_in.close()
hap1_auto_bam_in.close()
hap2_auto_bam_in.close()
hap1_bam_out.close()
hap2_bam_out.close()