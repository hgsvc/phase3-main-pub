from argparse import ArgumentParser 
import mhca
import os
from importlib import resources
import subprocess
import sys
import itertools
from pathlib import Path

def main(args):
    """ Main function extract full fasta from refseq file"""
    cwd = os.getcwd()

    # extract Refseq so you have the following fields
    ##bin	name	name2	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
    # example: 
    #802	NM_182701	GPX6	chr6	-	28503295	28515793	28504291	28515743	5	"28503295,28505702,28506311,28510750,28515656,"	"28504498,28505802,28506429,28510904,28515793,"	0	GPX6	cmpl	cmpl	"0,2,1,0,0,"

    with open(args.refseq_genes) as inf, open(args.outfile, "w") as out:
        for line in inf:
            if line.startswith('#'): continue

            _, accession, gene_name, chrom, strand, txstart, txend, cdsstart, cdsend, nrexons, exon_starts, exon_ends, _, _, cdsStartStat, cdsEndStat, exon_frames = line.rstrip().split()
            rc_string = ""
            if strand == "-": rc_string = "--reverse-complement"
            
            if cdsStartStat != "cmpl" or cdsEndStat != "cmpl": continue
            #if gene_name == "C4B_2": continue
            samtools_job = subprocess.run(["samtools", "faidx", args.reference, "chr6:"+ txstart + "-" + txend, "--length", "100000000", rc_string ], capture_output=True)
            output = samtools_job.stdout.decode("utf-8")
            #if len(output.split("\n")) > 3:
            #    print(output)
            idx, seq, _ = output.split("\n")[:3]
            exon_starts = [int(x) for x in exon_starts.strip("\"").rstrip(",").split(",")]
            exon_ends = [int(x) for x in exon_ends.strip("\"").rstrip(",").split(",")]
            #print(exon_starts)
            cdsend = int(cdsend)
            cdsstart = int(cdsstart)
            txstart = int(txstart)
            txend = int(txend)
            txlen = txend - txstart
            cdsstart_rel = cdsstart-txstart
            cdsend_rel = cdsend -txstart
            assert cdsstart >= txstart
            assert cdsend <= txend

            bounds_t = list(zip(exon_starts,exon_ends))
            bounds_rel = [(exstart-txstart, exstop-txstart) for exstart, exstop in bounds_t]
            nbounds = []
            start_found = False
            for b1, b2 in bounds_rel:
                if not start_found:
                    if cdsstart_rel < b2 and cdsstart_rel >= b1:
                        nbounds.append(b1)
                        nbounds.append(cdsstart_rel)
                        if cdsend_rel <=b2 and cdsend_rel > b1:
                            nbounds.append(cdsend_rel)
                            nbounds.append(b2)
                            break
                        else:
                            nbounds.append(b2)
                        start_found = True
                else:     
                    if cdsend_rel <=b2 and cdsend_rel > b1:
                        nbounds.append(b1)
                        nbounds.append(cdsend_rel)
                        nbounds.append(b2)
                        break
                    else:
                        nbounds.append(b1)
                        nbounds.append(b2)
            for b1, b2 in zip(nbounds[:-1],nbounds[1:]):
                if b1 > b2:
                    print(f"Problem with {gene_name}: {nbounds}")
                    break
            
            bounds = nbounds
            #bounds = [0] + bounds + [int(txend) - int(txstart)]
            #print(bounds)
            if strand == "-":
                bounds = [txlen-x-1 for x in bounds[::-1]]
            #print(bounds)
            full_fasta = ""
            segments = []
            for first, last in zip(bounds[0:-1], bounds[1:]):
                segments.append(seq[first+1:last+1])
                #full_fasta += seq[first:last] + "|"
            if segments[1][0:3] != "ATG": print(f"Not a proper start codon for {gene_name}: {segments[0]}")
            if segments[-2][-3:] not in {"TGA", "TAA", "TAG"}: print(f"Not a proper stop codon for {gene_name}: {segments[-2]}")
            full_fasta = "|".join(segments)
            out.write(">" + gene_name.replace("_","-") + "_tr1\n")
            out.write(full_fasta + "\n")

            

