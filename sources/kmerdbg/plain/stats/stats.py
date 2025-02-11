#! /usr/bin/env python

from __future__ import print_function
import argparse
import collections
import csv
import gzip
import sys
import datetime

# Parse command line
parser = argparse.ArgumentParser(description='K-mer summary')
parser.add_argument('-k', '--kmer', metavar='nonref.sample.gz', required=True, dest='kmer', help='k-mer file (required)')
parser.add_argument('-s', '--samples', metavar='samples.tsv', required=True, dest='sample', help='sample info (required)')
args = parser.parse_args()

pop = dict()
# Parse samples
if args.sample:
    with open(args.sample) as f:
        f_reader = csv.DictReader(f, delimiter="\t")
        for row in f_reader:
            pop[row['sample']] = row['pop']

# Count
singletons = collections.Counter()
doubletons = collections.Counter()
popspecific = collections.Counter()
totalkmer = collections.Counter()
num = 0
total = 0
if args.kmer:
    with gzip.open(args.kmer) as f:
        f_reader = csv.reader(f, delimiter="\t")
        for row in f_reader:
            if num > 1000000:
                total += 1
                now = datetime.datetime.now()
                print (now.strftime("%Y-%m-%d %H:%M:%S"), ": Processed", total, "million records", file=sys.stderr)
                num = 0
            num += 1
            samples = row[3].split(',')
            popsp = True
            if len(samples) == 1:
                singletons[samples[0]] += 1
            elif len(samples) == 2:
                if pop[samples[0]] == pop[samples[1]]:
                    doubletons[samples[0]] += 1
                    doubletons[samples[1]] += 1
                else:
                    popsp = False
            else:
                mypop = pop[samples[0]]
                for s in samples:
                    if mypop != pop[s]:
                        popsp = False
                        break
            for s in samples:
                totalkmer[s] += 1
                if popsp:
                    popspecific[s] += 1

print("sample", "pop", "kmer", "singleton", "doubleton", "popspecific", sep="\t")
for s in pop.keys():
    print(s, pop[s], totalkmer[s], singletons[s], doubletons[s], popspecific[s], sep="\t")
