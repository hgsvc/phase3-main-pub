# Inputs used for HGSVC3 analyses

## Pipeline

Config file used: config.yaml (stored in this folder)
 * merged_set: pav_variants_batch1_alt.vcf.gz, produced from PAV calls as described below
 * single_vcfs / callable regions: single-sample VCFs were donwloaded from: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/working/20231016_PAV_2.3.3/pav_2.3.3_hg38_mm2_vcf.tar.gz, callable regions from: https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/pav_callable_batch1.tar.gz
 * Illumina reads: "1000GP high-coverage cohorts" reads downloaded from EBI/ENA and merged into one file per sample
 * assembly_qvs: computed by Peter Ebert (stored in this folder)
 * regions: segmental duplications track for GRCh38 downloaded from UCSC
 * yak version: commit f37704a

## PAV data

merged PAV VCFs were downloaded from:

https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_snv_snv_alt.vcf.gz
https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_snv_snv_alt.vcf.gz.tbi

https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_indel_insdel_alt.vcf.gz
https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_indel_insdel_alt.vcf.gz.tbi

https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_sv_insdel_alt.vcf.gz
https://storage.googleapis.com/jax-beck-pub/hgsvc3/variant_calls/batch1/variants_batch1_sv_insdel_alt.vcf.gz.tbi

and merged into a single VCF using the command:

bcftools concat -a <files> | bgzip > pav_variants_batch1_alt.vcf.gz
