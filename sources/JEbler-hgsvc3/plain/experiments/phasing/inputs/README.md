Assembly stats computed and provided by Peter: hgsvc3_varcall_qvest.tsv


BED files obtained from:

GIAB stratifications (CHM13):
CHM13_notinalldifficultregions.bed: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/CHM13@all/Union/CHM13_notinalldifficultregions.bed.gz
CHM13_alllowmapandsegdupregions.bed: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/CHM13@all/Union/CHM13_alllowmapandsegdupregions.bed.gz
CHM13_AllTandemRepeatsandHomopolymers_slop5.bed: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/CHM13@all/LowComplexity/CHM13_AllTandemRepeatsandHomopolymers_slop5.bed.gz

CHM13 unique/shared regions:

chm13-unique.bed: Regions unique to CHM13 (not in GRCh38) Downloaded from UCSC "CHM13 unique" track for CHM13v2.0
chm13-shared-hg38.bed: Complement of unique regions, i.e. complement of the above BED. Generated with following commands:

cat chm13-unique.bed | sort -k1,1 -k2,2n > chm13-unique-sorted.bed
sort -k1,1 -k2,2n chm13v2.0.fa.fai > ref.fai
bedtools complement -i chm13-unique-sorted.bed -g ref.fai > chm13-shared-hg38.bed

