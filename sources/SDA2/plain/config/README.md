# SegDupAnnotation2 General Configuration

This workflow requires a `config/sd_analysis.json` configuration file. Please modify it to your needs. It is valudated against `workflow/schemas/config.schema.json` which also contains example parameters.

## Required Parameters

The key files required by this pipeline are an assembly file, a genemodel (for instance a refseq annotation file), and the PacBio reads bam.

| Parameter | Type | Description |
| --- | --- | ------ |
| species | String | The name, code, or identifier for the species or individual being processed. |
| reads | Array of Strings | List of paths to read files. Acceptable formats and extensions: bam, fastq.gz, fastq, fasta, and fa. (If none provided, workflow will run in 'bamless' mode and assume a read depth of 100 bp across the assembly. Output files will reflect this depth, and thus results may be confusing. For this reason it is not advised to use this workflow without pacbio read bams.) |
| asm | String | File path to assembly of this species/individual's genome in fasta format. If using haplotype-aware mode, place haplotype 1 here. |
| genemodel | String | Fasta filepath of chosen gene model. Assumes fasta headers have no space characters. |
| temp | String | File path to an existing directory for temporary files. |

## Optional Parameters

| Parameter | Type | Default | Description |
| --- | --- | --- | ------ |
| hap2 | String | N/A | File path to haplotype 2 of this species/individual's genome in fasta format. If provided, will activate haplotype-aware mode. The contig/scaffold/chromosome names in the hap2 file must be unique relative to the provided asm/hap1 file. |
| read_type | String | N/A | Note PacBio read technology type (CLR vs CCS). Currently used for metadata purposed only. |
| haploid_chrs | Array of Strings | N/A | List of haploid/sex chromosome names in given assembly. Will only be used if flag_autodetect_haloid_chrs is not explicitly set to 'true'. Only helpful in haplotype-unaware mode to ensure collapsed duplications on haploid chromosomes are detected appropriately. |
| flag_autodetect_haploid_chrs | Boolean | false | When true autodetects whether a chromosome is haploid or diploid to inform correct gene duplication accounting in rule F05_FindDups. Rule B05 detects haploid chromosomes by looking for chromosomes with half the mean depth of the assembly. Chromosomes smaller than 10 Mb will not be selected as haploid chromosomes. Only helpful in haplotype-unaware mode to ensure collapsed duplications on haploid chromosomes are detected appropriately. |
| bed_for_filtering_results | String | N/A | Bed file path of assembly coordinates to filter out of results. |
| flagger_components_to_filter | Array of Strings | N/A | Provide flagger bed file to bed_for_filtering_results option. Here list which flagger components to exclude from analysis. Reference https://github.com/mobinasri/flagger for component descriptions. |
| min_copy_identity | Number | 0.90 | Minimum gene copy identity to keep when the gene copy is compared to the original copy. |
| min_hit_length | Int | 5000 | Minimum hit length to keep in bases. |
| max_length_margin | Number | 0.10 | Keep gene copies with length within \<max_length_margin\> of the original gene's length. |
| min_gene_model_alignment | Number | 0.75 | Minimum percent alignment of gene model to gene copy as a decimal. |
| min_depth | Number | 0.05 | Minimum mean copy depth to keep as percentage of mean assembly depth. |
| flag_filt_single_exon_genes | Boolean | true | When true, keeps only genes with multiple exons. |
| flag_force_hmm | Boolean | false | When true, uses the hidden markov copy number caller instead of samtools' mpileup to determine depth in 100 bp reads. This is not a recommended option, and may fail if no reads map to a single chromosome/scaffold. |
| flag_assume_clear_and_unique_gene_codes | Boolean | true | When false, assumes gene model fasta headers are in default RefSeq format, and thus renames all headers based on gene symbol in parenthesis at end of header line. |
| flag_filt_uncharacterized_genes | Boolean | true | When true, filters out genes in gene model with gene names beginning with `LOC`. |
| uncharacterized_gene_name_prefix | String | N/A | Gene model gene name prefix to filter out when flag_filt_uncharacterized_genes is true. Also used by Network Filter to deprioritize uncharacterized genes when picking a representative gene for a gene family. |
| cluster_mem_mb_baby | Int | 1000 | The memory in MB a cluster node or cpu must provide for a computationally simple job. In practice this parameter is combined with a cluster_cpus_per_task_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_mem_mb_small | Int | 6000 | The memory in MB a cluster node or cpu must provide for a computationally simple job. In practice this parameter is combined with a cluster_cpus_per_task_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_mem_mb_medium | Int | 6000 | The memory in MB a cluster node or cpu must provide for a computationally mild job. In practice this parameter is combined with a cluster_cpus_per_task_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_mem_mb_large | Int | 6000 | The memory in MB a cluster node or cpu must provide for a computationally intense job. In practice this parameter is combined with a cluster_cpus_per_task_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_mem_mb_xlarge | Int | 6000 | The memory in MB a cluster node or cpu must provide for a computationally intense job. In practice this parameter is combined with a cluster_cpus_per_task_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_cpus_per_task_baby | Int | 1 | The number of cpus per task for a computationally simple rule. In practice this parameter is combined with a cluster_mem_mb_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_cpus_per_task_small | Int | 3 | The number of cpus per task for a computationally mild rule. In practice this parameter is combined with a cluster_mem_mb_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_cpus_per_task_medium | Int | 3 | The number of cpus per task for a computationally intense rule. In practice this parameter is combined with a cluster_mem_mb_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_cpus_per_task_large | Int | 3 | The number of cpus per task for a computationally intense rule. In practice this parameter is combined with a cluster_mem_mb_\<size\> parameter by some rules to create a SLURM or other cluster job. |
| cluster_runtime_short | Int | 240 | The walltime in minutes allocated for rules expected to take a relatively short amount of time (like 4 hrs). This parameter is only used if called by snakemake's --slurm command line paramter. |
| cluster_runtime_long | Int | 240 | The walltime in minutes allocated for rules expected to take a relatively long amount of time (like 24 hrs). This parameter is only used if called by snakemake's --slurm command line paramter. |
| override_mem | Int | -1 | Override the memory available in MB otherwise defined by the cluster_mem_mb_\<size\> parameters in MB. If set to -1, the cluster_mem_mb_\<size\> paramters will not be overwritten. |
| override_num_cores | Int | -1 | Override the number of allocated cores otherwise defined by the cluster_cpus_per_task_\<size\> parameters. If set to -1, the cluster_cpus_per_task_\<size\> paramters will not be overwritten. |
