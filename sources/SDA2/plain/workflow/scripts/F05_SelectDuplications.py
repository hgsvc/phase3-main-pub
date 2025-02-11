#!/usr/bin/env python

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 06/08/23

# Purpose: Pick out genes with resolved or collapsed duplications 
# (with copy num greater than 2 for chroms unless input as haploid/sex chromosome.)
# Input: resolved copies bed file path. Haploid chromosomes may also be specified.
# Output: Filtered input.

import sys
import argparse
import pathlib

# Parse Input
parser = argparse.ArgumentParser(prog='F05_SelectDuplications.py',
                                 description="Pick out genes with resolved or collapsed duplications.")
parser.add_argument("bed_filepath", type=pathlib.Path, help="Bed file with resolved copies and copy numbers. Must be pre-sorted by gene_name in column 4.")
parser.add_argument('-s','--sex_chrs_list', nargs='+', help="Space separated haploid/sex chromosomes list.")
parser.add_argument('--sex_chrs_list_filepath', type=pathlib.Path, help="Path of line separated haploid/sex chromosomes list. Supderseded by -s/--sex_chrs_list flag.")
parser.add_argument('--use_vcf_depth', action='store_true')
parser.add_argument('--phased_input', action='store_true', help="Assume input bed file has hits from a union of both haplotypes generated from a haplotype phased assembly.")
args=parser.parse_args()

# Vars for bed file index reference
CHR_I=0
START_I=1
END_I=2
GENE_I=3
OG_CHR_I=4
OG_START_I=5
OG_END_I=6
STRAND_I=7
HAPLOTYPE_I=8
P_IDENTITY_I=9
P_ACCURACY_I=10
COPY_I=11
EXON_SIZES_I=12
EXON_STARTS_I=13
GM_ALIGNMENT_I=14
DEPTH_RATIO_I=15
CN_BY_DEPTH_STDEV_I=16
CN_BY_DEPTH_I=17
CN_BY_VCF_I=18

if args.use_vcf_depth:
    CN_I=CN_BY_VCF_I
else:
    CN_I=CN_BY_DEPTH_I

EXPECTED_HAPLOID_CN=1
EXPECTED_DIPLOID_CN=2

# Create set for Sex Chr Lookup
if (args.sex_chrs_list):
    s=set(args.sex_chrs_list)
elif (args.sex_chrs_list_filepath):
    sexChrFile=open(args.sex_chrs_list_filepath)
    s=set()
    for line in sexChrFile:
        s.add(line.strip())
    sexChrFile.close()
else:
    s=set()

# Helper Function to print collapsed or multicopy genes
# Input non-empty 2D array of copies of a gene
def printDups(geneList):
    sexChrPresent=False
    geneCNCount=0
    hap1ResolvedCNCount=0
    hap2ResolvedCNCount=0
    for copy in geneList:
        if copy[CHR_I] in s:
            sexChrPresent=True
        if copy[HAPLOTYPE_I] == "haplotype1":
            hap1ResolvedCNCount+=1
        if copy[HAPLOTYPE_I] == "haplotype2":
            hap2ResolvedCNCount+=1
        geneCNCount+=int(copy[CN_I])
    if (sexChrPresent and geneCNCount > EXPECTED_HAPLOID_CN) or \
       (not sexChrPresent and geneCNCount > EXPECTED_DIPLOID_CN) or \
       (not args.phased_input and len(geneList) >= 2) or \
       (args.phased_input and len(geneList) >= 3) or \
       (args.phased_input and hap1ResolvedCNCount >= 2) or \
       (args.phased_input and hap2ResolvedCNCount >= 2):
        for copy in geneList:
            sys.stdout.write('\t'.join(copy)+'\n')

# Traverse bed file grouping genes
lastGene=""
geneList=[]
bedFile=open(args.bed_filepath)
for line in bedFile:
    copy=line.split()
    if copy[GENE_I]!=lastGene and geneList: # if genes not equal and geneList has elts
        printDups(geneList) # test if gene list can be printed if so print
        geneList=[]
    geneList.append(copy)
    lastGene=copy[GENE_I]
if geneList: # if geneList not empty
    printDups(geneList)
bedFile.close()
