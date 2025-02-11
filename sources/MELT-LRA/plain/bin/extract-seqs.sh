#!/bin/bash

# ----------------------------
# insertions
# ----------------------------

# INFO - extracted 15875/658478 sequences (~2.4%) from insertion variants
./extract_insertion_seqs.py --vcf=/local/projects-t3/DEVINElab/HGSVC3/pav_HG00514_CCS_SS_PG_PRR.vcf.gz --min_seqlen=50 > CCS-insertion-seqs-gt50.fsa

# INFO - extracted 17689/2025514 sequences (0.87%) from insertion variants
./extract_insertion_seqs.py --vcf=/local/projects-t3/DEVINElab/HGSVC3/pav_HG00514_CLR_SS_FLYE_PA.vcf.gz --min_seqlen=50 > CLR-insertion-seqs-gt50.fsa

# ----------------------------
# deletions
# ----------------------------

#INFO - extracted 10012/505660 sequences from deletion variants
#INFO - ignored 495634/505660 sequences below minimum length
#INFO - ignored 14/505660 sequences above maximum length
./extract_deletion_seqs.py --vcf=/local/projects-t3/DEVINElab/HGSVC3/pav_HG00514_CCS_SS_PG_PRR.vcf.gz --min_seqlen=50 --max_seqlen=50000 > CCS-deletion-seqs-gt50.fsa

#INFO - extracted 11408/673983 sequences from deletion variants
#INFO - ignored 662553/673983 sequences below minimum length
#INFO - ignored 22/673983 sequences above maximum length
./extract_deletion_seqs.py --vcf=/local/projects-t3/DEVINElab/HGSVC3/pav_HG00514_CLR_SS_FLYE_PA.vcf.gz --min_seqlen=50 --max_seqlen=50000 > CLR-deletion-seqs-gt50.fsa

