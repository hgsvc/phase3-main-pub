#!/usr/bin/env python

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 12/21/23

# Purpose: Take grouped isoforms table and merge groups that contain the same isoforms. Pick a consensus isoform based on greatest total depth, then length. Filter input bed file.
# Input: Grouped isoforms table.
# Output: Filtered bed file.

import sys
import argparse
import pathlib
import networkx as nx

# Parse Input
print("---- Running F06_mergeIsoFamsAndPickConsensus.py", file=sys.stderr)
print("-- Parsing Input", file=sys.stderr)
parser = argparse.ArgumentParser(description="Take grouped isoforms table and merge groups that contain the same isoforms.")
parser.add_argument("tsv_filepath", type=pathlib.Path, help="Input tab separated file of grouped isoforms. (The first column is a comma separated list of isoform codes without its location information. The second column has a comma separated list of isoform codes with each copy's location information.) (Only the first column is needed for this script.)")
parser.add_argument("tsv_out_filepath", type=pathlib.Path, help="Output tab separated file of the picked concensus and the isoforms in its group.")
parser.add_argument("bed_filepath", type=pathlib.Path, help="Input extended bed file produced by rule F05.")
args = parser.parse_args()

GENE_BED_I = 3
OG_START_BED_I = 5
OG_END_BED_I = 6
CN_BED_I = 17

# Create Graph
print("-- Creating Graph", file=sys.stderr)
# (nodes = gene name, edges = conections given in input tsv_filepath file)
gene_fam_graph = nx.Graph()

# Parse tsv Input File & Add To Graph
tsvFile=open(args.tsv_filepath)
for line in tsvFile:
    tsvLine=line.split()
    genes=tsvLine[0].split(',')
    for g in genes:
        gene_fam_graph.add_edge(genes[0],g)
tsvFile.close()

# Add node weights
bedFile=open(args.bed_filepath)
for line in bedFile:
    bedLine = line.split()
    gene_name = bedLine[GENE_BED_I]
    copy_num = 0
    if gene_fam_graph.has_node(gene_name):
        if 'copy_num' in gene_fam_graph.nodes[gene_name]:
            copy_num=gene_fam_graph.nodes[gene_name]['copy_num']
    copy_num += int(bedLine[CN_BED_I])
    gene_fam_graph.add_node(gene_name,copy_num=copy_num,og_length=int(bedLine[OG_END_BED_I])-int(bedLine[OG_START_BED_I]))
bedFile.close()

# Pick Consensus And Print
print("-- Picking Consensus", file=sys.stderr)

# Helper function for sorting gene labels (family groupings)
# by the number of total duplications and then its isoform's
# length.
# Accepts gene_name and returns (key,value) pair of (duplication
# count, original isoform length).
# Hardcoded to use gene_fam_graph
def get_dup_count_and_mean_len (gene_name):
    if 'copy_num' in gene_fam_graph.nodes[gene_name]:
        dup_count=gene_fam_graph.nodes[gene_name]['copy_num']
    else:
        dup_count=0
    if 'og_length' in gene_fam_graph.nodes[gene_name]:
        og_iso_len=gene_fam_graph.nodes[gene_name]['og_length']
    else:
        og_iso_len=0
    return (dup_count,og_iso_len)

coms_out=open(args.tsv_out_filepath,'w')
label_set=set()
for fam_nodes in nx.connected_components(gene_fam_graph):
    sub_nodes=sorted(list(fam_nodes),key=get_dup_count_and_mean_len, reverse=True)
    rep_name_picked=sub_nodes[0]
    label_set.add(rep_name_picked)
    print("\t".join([rep_name_picked,
                     ",".join(sub_nodes)]), file=coms_out)
coms_out.close()

# Filter Bed File
print("-- Filter bed file", file=sys.stderr)
bedIn=open(args.bed_filepath)
for line in bedIn:
    vals=line.split()
    if vals[GENE_BED_I] in label_set:
        print('\t'.join(vals), file=sys.stdout)
bedIn.close()