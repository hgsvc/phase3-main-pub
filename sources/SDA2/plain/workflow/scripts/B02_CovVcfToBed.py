#!/usr/bin/env python

import sys

inFile=open(sys.argv[1])
for line in inFile:
    if len(line) > 0 and line[0] == "#":
        continue
    vals=line.split()
    if vals[6] != "PASS":
        continue
    chrom=vals[0]
    start=vals[1]
    kvp=vals[7].split(";")
    end=kvp[1].split("=")[1]
    format=vals[-2].split(":")
    sample=vals[-1].split(":")
    #cn=sample[0]

    i=0
    cn="-"
    dp="-"
    for f in format:
        if f == "CN":
            cn=sample[i]
        if f == "DP":
            dp=sample[i]
        i+=1
    sys.stdout.write(chrom + "\t" + start + "\t" + end + "\t" + cn + "\t" + dp + "\n")
inFile.close()