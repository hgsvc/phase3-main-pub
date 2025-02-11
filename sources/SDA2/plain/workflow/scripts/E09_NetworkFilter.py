#!/usr/bin/env python

# Keon Rabbani
# krabbani@usc.edu
# Chaisson Lab
# 05/24/23

# Purpose: Identify consensus isoforms from isoform clusters using the Leiden community detection algorithm on hits whose exons overlap.
# Input: pafxe filepath, filepath to generate output communities tsv, prefix of uncharacterized gene names to deprioritize when picking conesnsus gene
# Output: First 8 cols of pafx input file with a final column noting 0 if filtered out or 1 if hit remains in final set.

### NOTE: At this time, input cannot be piped in, since input file is accessed multiple times.

import sys
import argparse
import pathlib
import numpy as np
import networkx as nx
import igraph as ig
import pandas as pd

# NETWORK TUNING PARAMS:
MIN_OVERLAP = 10 # remove edges in exon overlap graph if exon overlap is less than this parameter in bases

# Parse Input
print("---- Running E09_NetworkFilter.py", file=sys.stderr)
print("-- Parsing Input", file=sys.stderr)
parser = argparse.ArgumentParser(description="Identify consensus isoforms from isoform clusters using the Leiden community detection algorithm on hits whose exons overlap on the same strand.")
parser.add_argument("pafxe_filepath", type=pathlib.Path, help="Input pafxe file. (Including comma separated list of exon lengths, then a comma separated list of exon start positions.)")
parser.add_argument("communities_filepath", type=pathlib.Path, help="Output tab separated file of consensus isoform. (The first column is the consensus isoform. The second column has a comma separated list of isoforms in the consensus isoform's family.) The list in the 2nd column is sorted by largest to smallest isoforms in length.")
parser.add_argument("-u","--uncharacterized_gene_name_prefix", help="Prefix of uncharacterized genes. Used to deprioritize uncharacterized genes when identifying a consensus gene in a gene community/family.")
parser.add_argument("-f","--leiden_objective_function",choices=['CPM','modularity'],default="CPM",help="Objective function used by igraph's Leiden Community Detection Method implemented here. Enter either 'modularity' or 'CPM' for the Constant Potts Model. See igraph's community_leiden() documentation.")
parser.add_argument("-r","--leiden_resolution",default=0.25,type=float,help="The resolution parameter to use in the Leiden algorithm when forming gene family communities. Higher resolutions lead to more small communities, while lower resolutions lead to fewer larger communities. See igraph's community_leiden() documentation.")
parser.add_argument("-n","--leiden_n_iterations",default=10,type=int,help="The number of iterations to iterate the Leiden algorithm when forming gene family communities. Each iteration may improve the partition further. Using a negative number of iterations will run until a stable iteration is encountered. See igraph's community_leiden() documentation.")
args = parser.parse_args()

# Determine array dimensions
pafxeFile=open(args.pafxe_filepath)

GENE_PAFXE_I=0
STRAND_PAFXE_I=4
CHROM_PAFXE_I=5
START_PAFXE_I=7
END_PAFXE_I=8
EXON_SIZES_PAFXE_I=22
EXON_STARTS_PAFXE_I=23

numTotalCopies=0
numTotalExons=0
geneStrLen=0
chromStrLen=0
exonSizeStrLen=0
exonStartStrLen=0

for line in pafxeFile:
    pafxeLine=line.split()
    numTotalCopies+=1

    thisGeneStrLen=len(pafxeLine[GENE_PAFXE_I])
    if geneStrLen < thisGeneStrLen:
        geneStrLen = thisGeneStrLen
    
    thisChromStrLen=len(pafxeLine[CHROM_PAFXE_I])
    if chromStrLen < thisChromStrLen:
        chromStrLen = thisChromStrLen

    thisExonSizeStrLen=len(pafxeLine[EXON_SIZES_PAFXE_I])
    if exonSizeStrLen < thisExonSizeStrLen:
        exonSizeStrLen = thisExonSizeStrLen

    thisexonStartStrLen=len(pafxeLine[EXON_STARTS_PAFXE_I])
    if exonStartStrLen < thisexonStartStrLen:
        exonStartStrLen = thisexonStartStrLen

    numExons=len(pafxeLine[EXON_STARTS_PAFXE_I].split(","))
    numTotalExons+=numExons

pafxeFile.close()


# Create ndarrays
isoforms = np.empty([numTotalCopies],
                dtype={'names': ('gene', 'q_seq_len', 'q_start', 'q_end', 'strand', 'chrom', 't_seq_len',
                                'start', 'end', 'num_matching_bases_in_mapping', 'num_bases_in_mapping',
                                'mapping_qual', 'aln_type', 'num_chain_minimizers', 'chaining_score',
                                'chaining_score_of_secondary', 'seq_divergence',
                                'length_of_query_rgn_with_repetitive_seeds', 'haplotype', 'percent_identity',
                                'percent_accuracy','original','exon_sizes','exon_starts','gm_alignment',
                                'copy_id','picked','exon_sum'),
                        'formats':('U'+str(geneStrLen),'u8','u8','u8','U1','U'+str(chromStrLen),'i',
                                'u8','u8','u8','u8',
                                'u1','U1','U10','U10',
                                'U10','U10',
                                'U10','U11','d',
                                'd','U8','U'+str(exonSizeStrLen),'U'+str(exonStartStrLen),'d',
                                'u8','U3','u8')})

GENE_ISOFORMS_I=0
Q_SEQ_LEN_ISOFORMS_I=1
ORIGINAL_ISOFORMS_I=21
COPY_ID_ISOFORMS_I=25
PICKED_ISOFORMS_I=26
EXON_SUM_ISOFORMS_I=27

exons = np.empty([numTotalExons],
                 dtype={'names': ('copy_id', 'chrom', 'strand', 'start', 'end','sum'),
                        'formats': ('u8','U'+str(chromStrLen),'U1','u8','u8','u8')})

COPY_ID_EXONS_I=0
CHROM_EXONS_I=1
STRAND_EXONS_I=2
START_EXONS_I=3
END_EXONS_I=4
SUM_EXONS_I=5

# Fill ndarrays
pafxeFile=open(args.pafxe_filepath)

copyNum=0
exonNum=0
for line in pafxeFile:
    pafxeLine=line.split()
    exonSizes=[int(num) for num in pafxeLine[EXON_SIZES_PAFXE_I].split(",")]
    exonSum=sum(exonSizes)
    isoforms[copyNum]=(pafxeLine[GENE_PAFXE_I],int(pafxeLine[1]),int(pafxeLine[2]),int(pafxeLine[3]),pafxeLine[4],pafxeLine[CHROM_PAFXE_I],int(pafxeLine[6]),
                int(pafxeLine[START_PAFXE_I]),int(pafxeLine[END_PAFXE_I]),int(pafxeLine[9]),int(pafxeLine[10]),
                int(pafxeLine[11]),pafxeLine[12],pafxeLine[13],pafxeLine[14],
                pafxeLine[15],pafxeLine[16],
                pafxeLine[17],pafxeLine[18],float(pafxeLine[19]),float(pafxeLine[20]),pafxeLine[ORIGINAL_ISOFORMS_I],
                pafxeLine[EXON_SIZES_PAFXE_I],pafxeLine[EXON_STARTS_PAFXE_I],float(pafxeLine[24]),
                copyNum,'no',exonSum)
    exonStarts=[int(num) for num in pafxeLine[EXON_STARTS_PAFXE_I].split(",")]
    for i in range(len(exonSizes)):
        exonStart=int(pafxeLine[START_PAFXE_I])+exonStarts[i] # gene start plus exon start (exon start is relative to first position of gene, not chromosome)
        exonEnd=exonStart+exonSizes[i]
        exons[exonNum]=(copyNum,pafxeLine[CHROM_PAFXE_I],pafxeLine[STRAND_PAFXE_I],exonStart,exonEnd,exonSum) # copy_id, chrom, strand, start, end, sum
        exonNum+=1
    copyNum+=1

pafxeFile.close()

exons.sort(order=('chrom','start','end','copy_id'))

# Create Graph
print("-- Creating Graph", file=sys.stderr)
# (nodes = index of each gene, edges = exon overlap, weight overlap = num bases of overlap between the exons of two genes)
overlapGraph = nx.Graph()
# Add Nodes
for geneNum in range(numTotalCopies):
    overlapGraph.add_node(geneNum)

# Helper functions for adding edges.
# Function overlap: returns true if two genes overlap and false 
# otherwise. Assumes start position of first gene is less than
# start position of second gene.
# c1=chrom_gene1, r1=strand_gene1, e1=end_gene1,
# c2=chrom_gene2, r2=strand_gene2, s2=start_gene2
def overlap (c1,r1,e1,c2,r2,s2):
    return (c1==c2 and r1==r2 and s2<=e1)
    #return (c1==c2 and s2<=e1)
    # formerly this function was strand agnostic

# Function overlapDist: returns magnitude of overlap in bases of
# two genes. Assumes start position of first gene is less than
# start position of second gene and that both genes overlap.
# e1=end_gene1, s2=start_gene2, e2=end_gene2
def overlapDist (e1,s2,e2):
    if e1 < e2:
        return (e1-s2)
    else:
        return (e2-s2)

# Add Edges
i=0
for ex in exons:
    j=i+1
    while j<exons.size and overlap(ex[CHROM_EXONS_I],ex[STRAND_EXONS_I],ex[END_EXONS_I],exons[j][CHROM_EXONS_I],exons[j][STRAND_EXONS_I],exons[j][START_EXONS_I]):
        if overlapGraph.has_edge(exons[i][COPY_ID_EXONS_I],exons[j][COPY_ID_EXONS_I]):
            overlapSum=overlapGraph.get_edge_data(exons[i][COPY_ID_EXONS_I],exons[j][COPY_ID_EXONS_I])["overlap"] + overlapDist(ex[END_EXONS_I],exons[j][START_EXONS_I],exons[j][END_EXONS_I])
            overlapGraph.add_edge(exons[i][COPY_ID_EXONS_I],exons[j][COPY_ID_EXONS_I],overlap=overlapSum)
        elif ex[COPY_ID_EXONS_I] != exons[j][COPY_ID_EXONS_I]:
            overlapGraph.add_edge(exons[i][COPY_ID_EXONS_I],exons[j][COPY_ID_EXONS_I],overlap=overlapDist(ex[END_EXONS_I],exons[j][START_EXONS_I],exons[j][END_EXONS_I]))
        j+=1
    i+=1

# Calculating Heuristic Weights for each edge and pruning short overlaps
edgeList=list(overlapGraph.edges(data=True))
for u,v,o in edgeList:
    if o['overlap'] < MIN_OVERLAP:
        overlapGraph.remove_edge(u,v)
    else:
    # average of sum of exons of nodes connected by edge
        avg=((isoforms[u][EXON_SUM_ISOFORMS_I] + isoforms[v][EXON_SUM_ISOFORMS_I]) / 2)
        overlapOverAvg = o['overlap'] / avg
        sqOverlapOverAvg = overlapOverAvg * overlapOverAvg # square of average of sum of exons of nodes connected by edge
        overlapGraph.add_edge(u,v,overlap=o['overlap'],
                              overlapOverAvgLen=overlapOverAvg,
                              overlapOverAvgLen100=int(overlapOverAvg*100),
                              overlapOverAvgLen_min1=1-overlapOverAvg,
                              sqOverlapOverAvgLen=sqOverlapOverAvg)

del edgeList

print("-  "+str(len(overlapGraph.nodes))+" gene isoforms detected (nodes).", file=sys.stderr)
print("-  "+str(len(overlapGraph.edges))+" isoform to isoform exon overlaps (edges).", file=sys.stderr)

# Community Detection # Lieden
# https://igraph.org/python/doc/api/igraph.Graph.html#community_leiden
print("-- Detecting Communities", file=sys.stderr)
subgraphs=[]
communities_leiden=[]
for subgraph_nodes in nx.connected_components(overlapGraph):
    sub=overlapGraph.subgraph(subgraph_nodes).copy()
    subgraphs.append(sub)
    sub_ig = ig.Graph.from_networkx(sub)
    if len(list(sub.nodes)) > 1:
        cluster=sub_ig.community_leiden(weights="overlapOverAvgLen",objective_function=args.leiden_objective_function,resolution=args.leiden_resolution,n_iterations=args.leiden_n_iterations)
        for i,com in enumerate(cluster):
            com_nodes=[]
            for c in com:
                com_nodes.append(sub_ig.vs[c]['_nx_name'])
            communities_leiden.append(com_nodes)
    elif len(list(sub.nodes)) == 1:
        communities_leiden.append([list(sub.nodes)[0]])
print("-  "+str(len(communities_leiden))+" gene communities detected.", file=sys.stderr)

# Join Leiden Communities By Gene Names
# (nodes = gene name, edges = between communities_leiden family members, node weight ids = list of isoform ids)
gene_fam_graph = nx.Graph()

# Helper function for adding node to gene_fam_graph graph
# while maintaining node attribute. The node name is set
# to the gene_name. The first node attribute is a list of 
# isoform ids that all have the same node name. The second
# node attribute is the isoform_id of the 'Original' isoform
# copy.
# Returns gene_name
# Hardcoded to use isoforms ndarray and gene_fam_graph graph.
def add_node_with_attributes (isoform_id):
    gene_name = isoforms[isoform_id][GENE_ISOFORMS_I]
    if gene_fam_graph.has_node(gene_name): # update node
        gene_fam_graph.nodes[gene_name]['isoforms'].add(isoform_id)
    else: # create new node
        gene_fam_graph.add_node(gene_name, isoforms={isoform_id})
    # Record Originality
    if isoforms[isoform_id][ORIGINAL_ISOFORMS_I] == 'Original':
        gene_fam_graph.add_node(gene_name, original=isoform_id)
    return gene_name

for c in communities_leiden:
    last_node = next(iter(c))
    last_name = add_node_with_attributes(last_node)
    for n in c:
        name = add_node_with_attributes(n)
        gene_fam_graph.add_edge(last_name,name)
        last_node = n
        last_name = name

gene_fam_names = nx.connected_components(gene_fam_graph)

gene_families = []
for fam_nodes in gene_fam_names:
    subgraph = gene_fam_graph.subgraph(fam_nodes).copy()
    subgraph_ids=[]
    for n in list(subgraph.nodes()):
        #subgraph_ids.append(gene_fam_graph.nodes[n]['original']) # version to save original only
        subgraph_ids.extend(gene_fam_graph.nodes[n]['isoforms'])
    gene_families.append(subgraph_ids)
print("-  "+str(len(gene_families))+" gene families detected.", file=sys.stderr) 

# Pick Concensus Gene
print("-- Picking Concensus Isoform", file=sys.stderr)
# Helper function for deprioritizing uncharacterized genes.
# Function is_uncharacterized: returns true if prefix of
# uncharacterized genes provided and given gene name starts
# with given prefix.
def is_characterized (name):
    if not args.uncharacterized_gene_name_prefix:
        return True
    return not name.startswith(args.uncharacterized_gene_name_prefix)

# Helper function for sorting iso_names (family groupings)
# by isoform length.
# Accepts isoform_id and returns its original's length
def og_len (gene_name):
    og_id=next(iter(gene_fam_graph.nodes[gene_name]['isoforms']))
    return isoforms[og_id][Q_SEQ_LEN_ISOFORMS_I]

coms_out = open(args.communities_filepath,"w")
print("\t".join(["consensus_isoform","all_isoforms_in_family"]),file=coms_out)
for c in gene_families:
    # Identify Consensus Node
    big_n = next(iter(c))
    big_oa = 0
    big_n_characterized = False
    sub_c = overlapGraph.subgraph(c).copy()
    iso_names = []
    for n in c:
        sum_oa = 0
        n_characterized = is_characterized(isoforms[n][GENE_ISOFORMS_I])
        #if not ((not big_n_uncharacterized) and n_uncharacterized): # in the case where big_n is characterized and n is an uncharacterized gene, will always retain the old characterized gene as a possible consensus node allowing for the deprioritization of uncharacterized genes if an uncharacterized gene prefix is specified.
        if not (big_n_characterized and (not n_characterized)):
            for u,v in sub_c.edges(n): # for each edge of node n
                sum_oa += sub_c.get_edge_data(u,v)['overlapOverAvgLen']
            if (not big_n_characterized) and n_characterized: # prioritize characterized over uncharacterized genes
                big_oa = sum_oa
                big_n = n
                big_n_characterized = n_characterized
            else: # if both genes being compared have same characterization state prefer gene with larger oa
                if sum_oa > big_oa:
                    big_oa = sum_oa
                    big_n = n
                    big_n_characterized = n_characterized
        iso_names.append(isoforms[n][GENE_ISOFORMS_I])
    # Record Consensus Node
    picked_name=isoforms[big_n][GENE_ISOFORMS_I]
    for iso in gene_fam_graph.nodes[picked_name]['isoforms']:
        isoforms[iso][PICKED_ISOFORMS_I]='yes'
    iso_names=sorted(iso_names, key=og_len, reverse=True)
    print("\t".join([isoforms[big_n][GENE_ISOFORMS_I],
                     ",".join(iso_names)]),
          file=coms_out)
coms_out.close()

# Output Data
print("-- Printing Output", file=sys.stderr)
df=pd.DataFrame(isoforms)
df=df.drop('copy_id',axis=1) # delete copy_id col
df=df.drop('exon_sum',axis=1) # delete exon_sum col

df.to_csv(sys.stdout,sep="\t",encoding='utf-8',header=False,index=False)
