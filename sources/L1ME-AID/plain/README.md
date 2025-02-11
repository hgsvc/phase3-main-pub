# Current Version: L1ME-AID (v1.1.1 beta)
***If you are looking for an earlier version (e.g., v1.0.0-beta) check the previousVersions folder.*** </br>
-Note: Beginning at Version 1.1.0-beta L1ME-AID will check for TSD sequence if you provide the reference genome as -g and name your sequences as 'chromosome-position-anythingElse' (e.g., chr1-1002321-whatever). TSD check will only happen if you give -g a path to a reference file otherwise this functionality is skipped.

## L1ME-AID (L1 Mediated Annotation and Insertion Detector)
<p align="center">
<img src="limeaid.jpeg?raw=true" width="350" height="350">
</p>
L1ME-AID is a simple python script that takes two files as input:
  1) Fasta file of structural variation sequences you want to check for L1 mediated transposition events (Mobile element insertions <i>Alu</i> elements, L1, SVA, etc.)
  2) RepeatMasker .out file that was produced from running RepeatMasker on the Fasta file.

# Installation
Download the code and run the python file (see HOW TO RUN below).

# Requirements:
This program was tested on Linux (Mint 21.2). It has not been tested on Mac/PC.

## Hardware:
A standard computer with enough RAM to support the in-memory operations/dataset you provide.

## Python Dependencies:
  1) matplotlib
  2) pandas
  3) pysam
  4) os
  5) ast
  6) numpy
  7) json
  8) collections
  9) from tqdm import tqdm
  10) from Bio.Seq import Seq
  11) from Bio import SeqIO
  12) more_itertools 
  13) argparse

# HOW TO RUN:

## Minimum Code Necessary:
`python limeaid.py -i /your/path/to/fastafile/hgsvc_INS.fasta -r /your/path/to/repeatmaskerfile/hgsvc_INS.fasta.out -o /your/path/to/outfile/hgsvc_INS_limeaid.csv`
- Provide an input Fasta file (-i), the Repeatmasker .out file (-r), and the path to the csv that will be exported (-o).

## Demo:
Check out the demo file for a subset of 10 random calls from the previous HGSVC callset (demo.fasta), a repeatmasker file on more than the ten (1kgpCallset.fasta.out), and the L1ME-AID out file (demo.outFile):</br>

You can see if it works for you by running this code: `python limeaid.py -i demo.fasta -r demo.fasta.out -o ./demo.outFile`</br>
The expected run time is less than 5 seconds for this demo. You should expect to get the same file as the demo.outFile provided. 

## Other options

  a) -rd, Repeat Divergence, default=20% divergence; 20% maximum divergence for reading through all the repeatmasker annotations (I would leave this but you can change it by providing an integer (-rd 10)). <br>
  b) -d, divergences, default=6,15,15 ; 6% for <i>Alu</i> elements, 15% for L1s and 15% for SVAs (provide a each number in this list format x,y,z, EXAMPLE: -d 5,10,15 )<br>
  c) -l, length cutoffs, default=500,10000,10000; 500bp for <i>Alu</i> elements, 10kbp for L1, 10kbp for SVA (provide a each length in this list format x,y,z, EXAMPLE: -l 300,8000,3000) <br>
  d) -u, upper length cutoff, default=50000, if a sequence is longer than this don't even bother<br>
  e) -t, maximum tail start position, default=50; this will allow the start of an elements tail to begin at position 50 in the sequence (feel free to change as desired)<br>
  f) -yl, allowable LINE subfamilies, default=L1HS,L1PA2; These are the L1 subfamilies that L1ME-AID will consider to be active. If for example an L1PA3 is annotated this will be noted as OLDER LINE SUBFAMILY<br>
  g) -g, path to the reference genome utilized during SV calling (e.g., '/home/yourname/myfolder/hg38.fa'). *Note: Only provide a reference genome here if you want to have L1ME-AID try to call TSDs. If you do this you need to name the sequences in your fasta in a specific format (chromosome-position-anything else: chr1-101101-mysequence). L1ME-AID uses this chromosome and position information to look for the coordinate sequence that matches in your SV sequence. 

## Pay Attention
Please pay attention to the following columns:<br>
  a) ID = the IDs need to be the same in the fasta file as the repeatmasker out file<br>
  b) TE_Percentage: What proportion of the sequence was annotated as specific elements<br>
  c) TE_Designation: The type of element that makes up the largest proportion of the sequence<br>
  d) FILTER_RESULTS: This will provide you the results from the filters (divergences, lengths, etc.). It will say 'Good_Row' if the element passes the checks. HOWEVER, it will always say 'Good_Row' for sequences that are not Alu elements, L1s, and SVAs. Additionally, this DOES NOT take into account other columns (Tail_Type, Tail_Length, Unique_Element_Count).<br>
  e) Tail_Begins: Where in the sequence the tail begins (this number corresponds to the antisense orientation so we are numbering everything from left to right).<br>
  f) Tail_Type: This will try and tell you if there is a tail possibly what type is present. If you are looking for active mobile elements there should be a tail. FILTER_RESULTS will say 'Good_Row' even if there is NO TAIL. Make sure to filter your data by both columns.<br>
  g) Tail_Length: Tries to estimate the length of the tail. Finding tails uses a 5bp seed so the minimum should be 5bp.<br>
  h) Tail_Seed_Hits: How many unique 5bp seeds were found (not too important for the user).<br>
  i) Unique_Element_Count: This will tell you if multiple subfamilies were identified in a sequence (Ex: <i>AluS</i>, <i>AluY</i>). This is more important for <i>Alu</i> elements. The code doesn't care about this column but it's important for the user to decide. Usually it's not a big deal if L1s or SVAs have multiple subfamilies as RepeatMasker will sometimes annotate pieces of the element separately. <br>

# Please Remember:
This program will provide you all of the output no matter the results. You can always choose to ignore the results if you think something is wrongly annotated. It is only here to help you in finding active elements mobilized by L1 machinery. This is still in beta so it might not be perfect. I usually perform manual curation of calls after. 

# License
This project is covered under the GNU Lesser General Public License, version 3.0 (LGPL-3.0)

# Log
December 09, 2024: Fixed a bug where non-MEIs without TSDs were breaking the code. Version 1.1.1 was uploaded. Seems to work. </br>
November 19, 2024: Verson 1.1.0-beta uploaded, added the -g functionality to call TSDs that are present within the SV sequence and the reference genome. 
