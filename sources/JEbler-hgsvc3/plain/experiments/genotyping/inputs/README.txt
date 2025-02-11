
HPRC assembly gaps (CHM13):
CHM13_hprc-gaps.bed: Downloaded from: https://genome.cshlp.org/content/33/4/496/suppl/DC1, Supplementary Table S5



UCSC segmental duplications:

CHM13:
CHM13_SEDEF_segdups.bed: Downloaded from UCSC segmental duplications track (SEDEF) for CHM13v2.0



medically relevant genes (CMRG) analyzed in Wagner et al.:

CHM13_CMRG_wagner.bed: was created using these steps:

Step 1) 273 medically relevant genes studied in https://doi.org/10.1038/s41587-021-01158-1 were downloaded from: https://github.com/usnistgov/cmrg-benchmarkset-manuscript/blob/master/data/gene_coords/unsorted/GRCh38_mrg_bench_gene.bed and extracted with: cut -f4 GRCh38_mrg_bench_gene.bed > gene-names.tsv

Step 2) CHM13-based gene annotations were downloaded from: https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gff3.gz

Step 3) Coordinates in CHM13 space for the genes were extracted using the command: zcat Homo_sapiens-GCA_009914755.4-2022_07-genes.gff3 | python3 convert.py gene-names.tsv > CHM13_CMRG_wagner.bed



Other files:

CHM13_all.bed: full genome in CHM13 coordinates
CHM13_chrX.bed: full chrX chromosome in CHM13 coordinates
CHM13_chrY.bed: full chrY chromosome in CHM13 coordinates
