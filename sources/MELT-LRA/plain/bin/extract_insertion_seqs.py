#!/usr/bin/env python3

# Extract insertion sequences from PAV VCF and write them to FASTA.

import argparse
import csv
import gzip
import os
import re    
import sys

# increase CSV max field size
csv.field_size_limit(256 * 1024 * 1024)

def log_info(msg):
    sys.stderr.write("INFO - " + msg + "\n")
    sys.stderr.flush()

def log_fatal(msg):
    sys.stderr.write("FATAL - " + msg + "\n")
    sys.stderr.flush()
    sys.exit(1)

# ------------------------------------------------------
# extract_seqs
# ------------------------------------------------------

def extract_seqs(vcf, min_seqlen, max_seqlen):
    n_seqs = 0
    n_extracted = 0
    n_below_min = 0
    n_above_max = 0
    # check for exact duplicates
    vcf_ids = {}
    lnum = 0
    
    with gzip.open(vcf, 'rt') as fh:
        cr = csv.reader(fh, delimiter='\t')
        for row in cr:
            lnum += 1
            if re.match(r'^#', row[0]):
                continue
            
            (chrom, pos, vcf_id, ref, alt, qual, filt, info, fmt, *rest) = row

            inf_d = {}
            for inf in info.split(';'):
                (k, v) = inf.split('=')
                if k in inf_d:
                    log_fatal("key " + k + " already seen in " + info)
                inf_d[k] = v

            if inf_d['SVTYPE'] != 'INS':
                continue

            # strip reference base from alt sequence
            n_seqs += 1
            
            if len(ref) != 1:
                log_fatal("len(REF) != 1")
            if alt[0] != ref:
                log_fatal("alt[0] (" + alt[0] + " != ref (" + ref + ")")

            alt = alt[1:]
            if len(alt) >= min_seqlen and len(alt) <= max_seqlen:
                n_extracted += 1
                if vcf_id in vcf_ids:
                    log_fatal("duplicate vcf_id " + vcf_id + " at line " + str(lnum))
                vcf_ids[vcf_id] = True
                print(">", vcf_id + " " + chrom + ":" + pos + " " + info)
                print(alt)
            elif len(alt) < min_seqlen:
                n_below_min += 1
            elif len(alt) > max_seqlen:
                n_above_max += 1

        log_info("extracted " + str(n_extracted) + "/" + str(n_seqs) + " sequences from insertion variants")
        log_info("ignored " + str(n_below_min) + "/" + str(n_seqs) + " sequences below minimum length")
        log_info("ignored " + str(n_above_max) + "/" + str(n_seqs) + " sequences above maximum length")
    
# ------------------------------------------------------
# main()
# ------------------------------------------------------

def main():

    # input
    parser = argparse.ArgumentParser(description='Extract insertion sequences.')
    parser.add_argument('--vcf', required=True, help='Path to input VCF file.')
    parser.add_argument('--min_seqlen', required=False, type=int, default=1, help='Minimum insertion sequence length.')
    parser.add_argument('--max_seqlen', required=False, type=int, default=1000000, help='Maximum insertion sequence length.')
    args = parser.parse_args()
    extract_seqs(args.vcf, args.min_seqlen, args.max_seqlen)
        
if __name__ == '__main__':
    main()



