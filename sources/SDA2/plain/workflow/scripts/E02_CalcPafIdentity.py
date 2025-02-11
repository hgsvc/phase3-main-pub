#!/usr/bin/env python

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 05/11/23

# Purpose: Determine paf hit sequence identities through Needleman-Wunsch alignment.
# Input: sequence assembly file path, temp dir path, and paf file path (contains exactly 1 line) in this order
# Output format has the following 21 columns: query_seq_name, query_seq_length, query_start, query_end, strand,
#    target_seq_name, target_seq_length, target_start, target_end, num_matching_bases_in_mapping, num_bases_in_mapping,
#    mapping_quality, alignment_type, num_chain_minimizers, chaining_score, chaining_score_of_secondary, seq_divergence,
#    length_of_query_rgn_with_repetitive_seeds, percent_identity, percent_accuracy, edlib_cigar
#    -> refer to minimap2 documentation for details on first 18 columns

import sys
import argparse
from Bio import SeqIO
import edlib
import pysam

# Parse Input
parser = argparse.ArgumentParser(description="Determine paf hit sequence identities through Needleman-Wunsch alignment.")
#parser.add_argument("pafx_filepath", help="Input pafx file. Pafx is first 12 columns from paf format followed by 8 columns defined in script E02_CalcPafIdentity.py. Then followed by line number of pafx file, then a zero.")
parser.add_argument("asm_filepath", help="Fasta file containing query and target sequences.")
parser.add_argument("tmp_dirpath", help="Temporary directory filepath.")
parser.add_argument("pafLine_dirpath", help="paf file path (contains exactly 1 line) containing query and target coordinates.")
args = parser.parse_args()

# Global Variables
Q_NAME_I=0
Q_LEN_I=1
Q_START_I=2
Q_END_I=3
STRAND_I=4
T_NAME_I=5
T_LEN_I=6
T_START_I=7
T_END_I=8
NUM_MATCHING_I=9
MAP_LEN_I=10
MAP_QUAL_I=11
PAF_TAGS_I=12

# Extract Paf Line
pafIn=open(args.pafLine_dirpath)
pafLine=pafIn.readline()
pafIn.close()

# Parse Paf Line
paf_vals=pafLine.split(maxsplit=12)
tp='-'
cm='-'
s1='-'
s2='-'
dv='-'
rl='-'
hp='-'
if 0 <= PAF_TAGS_I < len(paf_vals):
    for val in paf_vals[PAF_TAGS_I].split():
        if val.startswith('tp'):
            tp=val.split(sep=':')[2]
        elif val.startswith('cm'):
            cm=val.split(sep=':')[2]
        elif val.startswith('s1'):
            s1=val.split(sep=':')[2]
        elif val.startswith('s2'):
            s2=val.split(sep=':')[2]
        elif val.startswith('dv'):
            dv=val.split(sep=':')[2]
        elif val.startswith('rl'):
            rl=val.split(sep=':')[2]
        elif val.startswith('hp'):
            hp=val.split(sep=':')[2]

# Cut Seqs
srcRgn=paf_vals[Q_NAME_I].split(sep='/')[1]
pysam.faidx(args.asm_filepath,srcRgn,'-o',args.tmp_dirpath+'/src_rgn.fasta')
srcSeq=SeqIO.read(args.tmp_dirpath+'/src_rgn.fasta','fasta').seq

trgRgn=paf_vals[T_NAME_I]+':'+paf_vals[T_START_I]+'-'+paf_vals[T_END_I]
pysam.faidx(args.asm_filepath,trgRgn,'-o',args.tmp_dirpath+'/trg_rgn.fasta')
trgSeq=SeqIO.read(args.tmp_dirpath+'/trg_rgn.fasta','fasta').seq
if paf_vals[STRAND_I]=="-":
    trgSeq=trgSeq.reverse_complement()

# Alignment
alignment = edlib.align(srcSeq,trgSeq,mode="NW",task="path")

# Calculate Stats & Print
if alignment.get('cigar') is None:
    # if no alignment, set null stats
    pIdentity=0
    pAccuracy=0
    empty_cigar=str(alignment.get('editDistance'))+"X"

    # Print Identity & Accuracy Nulls with paf
    sys.stdout.write('\t'.join(paf_vals[:12])+'\t'+tp+'\t'+str(cm)+'\t'+str(s1)+'\t'+str(s2)+'\t'+str(dv)+'\t'+str(rl)+'\t'+str(hp)+'\t'+str(pIdentity)+'\t'+str(pAccuracy)+'\t'+empty_cigar+'\n')
else:
    # Count matches, mismatches, and gaps
    nMatch=0
    nMisMatch=0
    nInDelEvent=0
    nInDel=0
    lastC=""
    lastNum=""

    for c in alignment.get('cigar'):
        if c.isdigit():
            lastNum+=c
        else:
            if (c=='='):
                nMatch+=int(lastNum)
            elif (c=='X'):
                nMisMatch+=int(lastNum)
            elif (c=='I' or c=='D'):
                if (lastC!='I' and lastC!='D'):
                    nInDelEvent+=1
                nInDel+=int(lastNum)
            lastNum=""
            lastC=c
            # Note: edlib.align only uses the following 4 symbols in the cigar string: '=IDX'

    # Print Identity & Accuracy with paf
    pIdentity=nMatch/(nMatch+nMisMatch+nInDelEvent)
    pAccuracy=nMatch/(nMatch+nMisMatch+nInDel)

    sys.stdout.write('\t'.join(paf_vals[:12])+'\t'+tp+'\t'+str(cm)+'\t'+str(s1)+'\t'+str(s2)+'\t'+str(dv)+'\t'+str(rl)+'\t'+str(hp)+'\t'+str(pIdentity)+'\t'+str(pAccuracy)+'\t'+alignment.get('cigar')+'\n')