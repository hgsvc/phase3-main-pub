# PAV merged set analysis

Pipeline that computed QV values from merged / single sample callsets. 

## Inputs

The pipeline expects the following inputs:
* **merged callset VCF**: path to a multi-sample VCF containing the merged variant calls
* **single-sample VCFs**: one callset VCF for each sample/haplotype (calls before merging)
* **callable regions**: BED file for each haplotype defining callable regions
* **Illumina reads**: Illumina reads for all samples
* **reference genome**: reference genome FASTA
* **Optional: precomputed assembly QVs**

The location of these files must be specified in the `` config/config.yaml `` file. It contains the following fields:

```bat

# path to merged, multi-sample VCF
merged_set: "path/to/merged.vcf.gz"

# TSV file specifying paths to single-sample VCFs
# format:   <sample>   </path/to/calls.vcf.gz>  </path/to/callable-regions.bed> 
single_vcfs: "/path/to/single-sample-files.tsv"

# TSV file specifying paths to FASTA/Q reads
# format: <sample>   </path/to/reads.fa>
reads: "/path/to/reads.tsv"

# reference genome underlying calls
reference: "/path/to/reference.fa"

# OPTIONAL: precomputed QVs for assemblies
# format: <H1|H2> <sample-name> <QV value>.
# If not available, set to: ""
assembly_qvs: "/path/to/assembly_qvs.tsv"

# path to yak executable
yak: "./yak/yak"

```

## Outputs

The QV values and plots will be written to: `` results/evaluation/ ``


## What the pipeline does

For each haplotype, the pipeline implants variants from the merged and single-sample VCFs into the reference genome such that each interval of callable regions becomes a contig. From these contigs, QV values are estimated.


## How to run the pipeline

Prepare the config file as explained above. Then, run the pipeline using the command:

``` bat
snakemake --use-conda -j <number of cores>
```
