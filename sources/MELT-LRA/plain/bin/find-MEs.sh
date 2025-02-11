#!/bin/tcsh

# run on everything except chrY
./find_MEs.py --vcf=./HGSVC3/pav_HG00514_CCS_SS_PG_PRR.vcf.gz \
 --fasta_dir=./refs/hg38.analysisSet.chroms \
 --alu_water=CCS-ALU-water.out.gz \
 --alu_water_rev=CCS-ALU-water-rev.out.gz \
 --sva_water=CCS-SVA-water.out.gz \
 --sva_water_rev=CCS-SVA-water-rev.out.gz \
 --line_water=CCS-LINE1-water.out.gz \
 --line_water_rev=CCS-LINE1-water-rev.out.gz \
 --min_seqlen=100 --min_pctid=90 --min_pctcov=85 --skip_seqids='chrY'  >all-MEs.txt 
exit

# run on chr1 insertions only
# ~20s with chr1 only VCF
#gzcat ./HGSVC3/pav_HG00514_CCS_SS_PG_PRR.vcf.gz | perl -ne 'print if (/^(##|chr1\t)/);' > pav-chr1-only.vcf
#gzip pav-chr1-only.vcf
./find_MEs.py --vcf=./pav-chr1-only.vcf.gz \
 --fasta_dir=./refs/hg38.analysisSet.chroms \
 --alu_water=CCS-ALU-water.out.gz \
 --alu_water_rev=CCS-ALU-water-rev.out.gz \
 --sva_water=CCS-SVA-water.out.gz \
 --sva_water_rev=CCS-SVA-water-rev.out.gz \
 --line_water=CCS-LINE1-water.out.gz \
 --line_water_rev=CCS-LINE1-water-rev.out.gz \
 --min_seqlen=100 --min_pctid=90 --min_pctcov=85 --seqid=chr1
exit


