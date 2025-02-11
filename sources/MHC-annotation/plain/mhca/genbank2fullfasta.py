from lxml import etree
from argparse import ArgumentParser 
from Bio import GenBank, SeqIO
from Bio.Seq import Seq
from collections import Counter, defaultdict
import sys
import re
from importlib import resources
import mhca.data


parser = ArgumentParser() 
parser.add_argument("input_gb")
parser.add_argument("output_s2s_fa")
parser.add_argument("output_full_fa")
args = parser.parse_args()
transcripts = Counter()
seen_transcripts = defaultdict(dict)
gene2locus = {}

with resources.open_text(mhca.data, 'LRG_RefSeqGene.tsv') as inf:
    for line in inf: 
        if line.startswith('#'): continue
        _, _, gene, version= line.rstrip().split()[:4]
        gene2locus[gene] =version 

with open(args.input_gb) as inf, open(args.output_s2s_fa, "w") as outs2s, open(args.output_full_fa, "w") as outfull:
    for record in GenBank.parse(inf):
                    
        print(record.version)
        gene = "noname"
        for feat in record.features:
            if feat.key == "CDS":
                #print(feat.location)
                for qual in feat.qualifiers:
                    if qual.key == "/gene=":
                        gene = qual.value.strip("\"")
                if gene in gene2locus and gene2locus[gene] != record.version: 
                    print(f"skipped gene {gene} in locus {record.version}")
                    continue

                if ">" in feat.location or "<" in feat.location: continue # don't want fuzzy locations
                if feat.location in seen_transcripts[gene]: 
                    print(f"skipped transcript for gene {gene}")
                    continue
                else: seen_transcripts[gene][feat.location] = record.sequence

                m= re.findall("(\d+..\d+)", feat.location)
                start = int(m[0].split("..")[0])
                end = int(m[-1].split("..")[1])
                transcripts[gene] += 1
                ct = transcripts[gene]

                outfull.write(">" + gene + "_tr" + str(ct) + "\n")
                seq = record.sequence
                lastend = 0
                temp_seq = ""
                for ion in m:
                    exons,exone = [int(x) for x in ion.split("..")]
                    temp_seq += seq[lastend:exons-1] + "|"
                    temp_seq += seq[exons-1:exone] + "|"
                    lastend = exone
                temp_seq += seq[lastend:]
                if feat.location.startswith("complement"):
                    out_ar = [str(Seq(x).reverse_complement()) for x in temp_seq.split("|")]
                    out_ar.reverse()
                    out_seq = "|".join(out_ar)
                else:
                    out_seq = temp_seq
                outfull.write(out_seq + "\n")

    print("legnth of seen_tr: " + str(len(seen_transcripts)))
    for gene,transcripts in seen_transcripts.items():
        smallest = -1
        biggest = 0
        print(gene)
        for transcript, sequence in transcripts.items():
            m= re.findall("(\d+..\d+)", transcript)
            start = int(m[0].split("..")[0])
            end = int(m[-1].split("..")[1])
            if start < smallest or smallest == -1: smallest = start
            if end > biggest: biggest= end
        if transcript.startswith("complement"):
            out_seq = str(Seq(sequence)[smallest-1:biggest].reverse_complement())
        else:
            out_seq = str(Seq(sequence)[smallest-1:biggest])
        outs2s.write(">" + gene + "\n")
        outs2s.write("".join(out_seq) + "\n")

#print(seen_transcripts.keys())
#for tr in seen_transcripts["LST1"]:
#    print(tr)
