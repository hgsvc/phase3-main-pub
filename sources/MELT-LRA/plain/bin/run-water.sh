#!/bin/bash

export WATER_ARGS='-gapopen 10.0 -gapextend 0.5'

# ALU - 281 bp  (253 =~ 90%)
# SVA - 1316 bp (1184 =~ 90%)
# LINE1 - 6019 bp (5417 ~= 90%)

# ----------------------------
# insertions
# ----------------------------

# CCS
water -asequence ALU.fa -bsequence CCS-insertion-seqs-gt50.fsa $WATER_ARGS -outfile CCS-ALU-water.out
water -asequence LINE1.fa -bsequence CCS-insertion-seqs-gt50.fsa $WATER_ARGS -outfile CCS-LINE1-water.out
water -asequence SVA.fa -bsequence CCS-insertion-seqs-gt50.fsa $WATER_ARGS -outfile CCS-SVA-water.out

# CLR
water -asequence ALU.fa -bsequence CLR-insertion-seqs-gt50.fsa $WATER_ARGS -outfile CLR-ALU-water.out
water -asequence LINE1.fa -bsequence CLR-insertion-seqs-gt50.fsa $WATER_ARGS -outfile CLR-LINE1-water.out
water -asequence SVA.fa -bsequence CLR-insertion-seqs-gt50.fsa $WATER_ARGS -outfile CLR-SVA-water.out

# ----------------------------
# deletions
# ----------------------------

# CCS
water -asequence ALU.fa -bsequence CCS-deletion-seqs-gt50.fsa $WATER_ARGS -outfile CCS-ALU-water-del.out
water -asequence LINE1.fa -bsequence CCS-deletion-seqs-gt50.fsa $WATER_ARGS -outfile CCS-LINE1-water-del.out
water -asequence SVA.fa -bsequence CCS-deletion-seqs-gt50.fsa $WATER_ARGS -outfile CCS-SVA-water-del.out

# CLR
water -asequence ALU.fa -bsequence CLR-deletion-seqs-gt50.fsa $WATER_ARGS -outfile CLR-ALU-water-del.out
water -asequence LINE1.fa -bsequence CLR-deletion-seqs-gt50.fsa $WATER_ARGS -outfile CLR-LINE1-water-del.out
water -asequence SVA.fa -bsequence CLR-deletion-seqs-gt50.fsa $WATER_ARGS -outfile CLR-SVA-water-del.out
