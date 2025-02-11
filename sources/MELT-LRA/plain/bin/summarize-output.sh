#!/bin/tcsh

# ----------------------------
# insertions
# ----------------------------

## insertions (except for 14/505660 above 50kb)

#n_lines=674492
#n_alignments=17689
#n >= 90% identity and >= 90% coverage = 630
./summarize_water_output.py --fasta=CLR-insertion-seqs-gt50.fsa --water=CLR-ALU-water.out.gz --min_pctid=90 --min_pctcov=90
# >=80% id >=0% cov = 1320
# >=85% id >=0% cov = 1107
# >=87% id >=0% cov = 971
# >=90% id >=0% cov = 850

# >=85% id >=85% cov = 792

#n_lines=1674350
#n_alignments=17689
#n >= 90% identity and >= 90% coverage = 25
./summarize_water_output.py --fasta=CLR-insertion-seqs-gt50.fsa --water=CLR-LINE1-water.out.gz --min_pctid=90 --min_pctcov=90
#n_lines=1180154
#n_alignments=17689
#n >= 90% identity and >= 90% coverage = 0
# longest is 567
./summarize_water_output.py --fasta=CLR-insertion-seqs-gt50.fsa --water=CLR-SVA-water.out.gz --min_pctid=90 --min_pctcov=90

# ----------------------------
# deletions
# ----------------------------

## deletions (except for 22/673983 above 50kb)
#n_lines=440131
#n_alignments=11408
#n >= 90% identity and >= 90% coverage = 561
./summarize_water_output.py --fasta=CLR-deletion-seqs-gt50.fsa --water=CLR-ALU-water-del.out.gz --min_pctid=90 --min_pctcov=90
#n_lines=1008939
#n_alignments=11408
#n >= 90% identity and >= 90% coverage = 27
./summarize_water_output.py --fasta=CLR-deletion-seqs-gt50.fsa --water=CLR-LINE1-water-del.out.gz --min_pctid=90 --min_pctcov=90
#n_lines=715067
#n_alignments=11408
#n >= 90% identity and >= 90% coverage = 0
./summarize_water_output.py --fasta=CLR-deletion-seqs-gt50.fsa --water=CLR-SVA-water-del.out.gz --min_pctid=90 --min_pctcov=90
