#!/usr/bin/env python

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 02/20/24

# Purpose: Phase aligned reads. 
# Assumption: input bams are in lexicographic order, and PG tag already exists in input bam headers
# Input: Path to bam of reads aligned to haplotype1 and path to bam of reads aligned to haplotype2.
# Output: 

import sys
import argparse
import pathlib
import pysam

### Parse Input
print("---- Running A06H_partitionReads.py", file=sys.stderr)
print("-- Parsing Input", file=sys.stderr)
parser = argparse.ArgumentParser(description="Partition reads given their mappings to each haplotype.")
parser.add_argument("hap1_bam", type=pathlib.Path, help="Path to bam of reads aligned to hap1.")
parser.add_argument("hap2_bam", type=pathlib.Path, help="Path to bam of reads aligned to hap2.")
parser.add_argument("hap1_only_bam", type=pathlib.Path, help="Output bam filepath for reads phased to hap1.")
parser.add_argument("hap2_only_bam", type=pathlib.Path, help="Output bam filepath for reads phased to hap2.")
parser.add_argument("hap1_auto_bam", type=pathlib.Path, help="Output bam filepath for perfectly autozygous reads aligned to hap1.")
parser.add_argument("hap2_auto_bam", type=pathlib.Path, help="Output bam filepath for perfectly autozygous reads aligned to hap2.")
parser.add_argument("unmapped_bam", type=pathlib.Path, help="Output bam filepath for unmapped reads (technically aligned to hap1).")
args = parser.parse_args()

### Open Files and Generate Headers
print("-- Open Files and Generate Headers", file=sys.stderr)
## Open Input Files
hap1_bam_in = pysam.AlignmentFile(args.hap1_bam, "rb")
hap2_bam_in = pysam.AlignmentFile(args.hap2_bam, "rb")

## Create Header Dictionaries
hap1_header_dict_updated = hap1_bam_in.header.as_dict()
hap2_header_dict_updated = hap2_bam_in.header.as_dict()
hap1_bam_last_PG_ID = hap1_header_dict_updated["PG"][-1]["ID"]
hap2_bam_last_PG_ID = hap2_header_dict_updated["PG"][-1]["ID"]
hap1_new_PG_entry = {'ID':'A06H_partitionReads', 'PN':'A06H_partitionReads.py', 'PP':hap1_bam_last_PG_ID,'CL':" ".join(sys.argv)}
hap2_new_PG_entry = {'ID':'A06H_partitionReads', 'PN':'A06H_partitionReads.py', 'PP':hap2_bam_last_PG_ID,'CL':" ".join(sys.argv)}
hap1_header_dict_updated["PG"].append(hap1_new_PG_entry)
hap2_header_dict_updated["PG"].append(hap2_new_PG_entry)

## Open Output Files
hap1_bam_out      = pysam.AlignmentFile(args.hap1_only_bam, 'wb', header=hap1_header_dict_updated)
hap2_bam_out      = pysam.AlignmentFile(args.hap2_only_bam, 'wb', header=hap2_header_dict_updated)
hap1_auto_bam_out = pysam.AlignmentFile(args.hap1_auto_bam, 'wb', header=hap1_header_dict_updated)
hap2_auto_bam_out = pysam.AlignmentFile(args.hap2_auto_bam, 'wb', header=hap2_header_dict_updated)
umap_bam_out      = pysam.AlignmentFile(args.unmapped_bam,  'wb', header=hap1_header_dict_updated)

### Iterate through bams
print("-- Iterate through bams", file=sys.stderr)
hap1_iter = hap1_bam_in.fetch(until_eof=True)
hap2_iter = hap2_bam_in.fetch(until_eof=True)

# Helper function for iterating through reads.
# Accepts pysam iterator and returns next primary or unmapped read.
# Thus secondary and supplemental reads are skipped. # TODO
def get_next_read(hap_iter):
    next_read = next(hap_iter, None)
    # while read exists, is not primary, and is not unmapped
    while next_read and not (next_read.flag & 0x900 == 0) and not next_read.is_unmapped:
        next_read = next(hap_iter, None)
    return next_read

hap1_read = get_next_read(hap1_iter)
hap2_read = get_next_read(hap2_iter)

# Assumes bam file is in lexicographic order
while hap1_read and hap2_read:
    if hap1_read.query_name == hap2_read.query_name:

        if hap1_read.is_unmapped and hap2_read.is_unmapped: # Segregate mutually unmapped reads
            umap_bam_out.write(hap1_read)
        elif hap1_read.is_mapped and hap2_read.is_mapped: # Pick 
            hap1_read_score = hap1_read.get_tag('AS')
            hap2_read_score = hap2_read.get_tag('AS')
            if (hap1_read_score > hap2_read_score):
                hap1_bam_out.write(hap1_read)
            elif (hap1_read_score < hap2_read_score):
                hap2_bam_out.write(hap2_read)
            else: # if DP Alignment Scores are exactly equal
                hap1_auto_bam_out.write(hap1_read)
                hap2_auto_bam_out.write(hap2_read)
        elif hap1_read.is_mapped:
            hap1_bam_out.write(hap1_read)
        else: #hap2_read.is_mapped
            hap2_bam_out.write(hap2_read)

        hap1_read = get_next_read(hap1_iter)
        hap2_read = get_next_read(hap2_iter)
    elif hap1_read.query_name < hap2_read.query_name:
        if hap1_read.is_unmapped:
            umap_bam_out.write(hap1_read)
        else:
            hap1_bam_out.write(hap1_read)
        hap1_read = get_next_read(hap1_iter)
    elif hap2_read.query_name < hap1_read.query_name:
        if hap2_read.is_unmapped:
            umap_bam_out.write(hap2_read)
        else:
            hap2_bam_out.write(hap2_read)
        hap2_read = get_next_read(hap2_iter)
    
while hap1_read:
    if hap1_read.is_unmapped:
        umap_bam_out.write(hap1_read)
    else:
        hap1_bam_out.write(hap1_read)
    hap1_read = get_next_read(hap1_iter)

while hap2_read:
    if hap2_read.is_unmapped:
        umap_bam_out.write(hap2_read)
    else:
        hap2_bam_out.write(hap2_read)
    hap2_read = get_next_read(hap2_iter)

### Close Files
print("-- Close Files", file=sys.stderr)
hap1_bam_in.close()
hap2_bam_in.close()
hap1_bam_out.close()
hap2_bam_out.close()
hap1_auto_bam_out.close()
hap2_auto_bam_out.close()
umap_bam_out.close()