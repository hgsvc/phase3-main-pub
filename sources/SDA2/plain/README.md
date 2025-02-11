# Snakemake Workflow: `SegDupAnnotation2`

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.25.0-brightgreen.svg)](https://snakemake.github.io)

## Overview

A snakemake workflow for counting gene duplications given PacBio reads, a genome assembly, and a gene model.  
Successor to [SegDupAnnotation](https://github.com/ChaissonLab/SegDupAnnotation).

:exclamation: Caution: This workflow is under active development. The main branch should run without any errors, but feel free reach out for assistance or submit a bug report in case you run into any problems or notice any errors.

## Usage
Usage guide is listed in order of ease of use and portability.  
Internet access is required for snakemake to download the docker image and conda environments.

### Singularity
- Install Snakemake and Singularity.
- Clone github repository.
- Create configuration file per config/README.md specifications.
    - If necessary don't forget to use the 'override_mem' and 'override_num_cores' parameters.
- From repository root, run (which will automatically download and run the dockerfile from within singularity):
    - `snakemake -c 1 -j 250 --use-singularity --use-conda --singularity-args " --bind \<path to an input file\>\[,\path to another input file\] "`
    - Don't forget to bind paths of input files per the given config file.

### Bare Metal with Conda
- Install Snakemake and Mamba.
- Clone github repository.
- Create configuration file per config/README.md specifications.
- From repository root, run (which will automatically download and use pre-defined conda environments):
    - `snakemake -c 1 -j 250 --use-conda`

### Bare Metal with Conda and SLURM
- Install Snakemake and Mamba.
- Clone github repository.
- Create configuration file per config/README.md specifications.
- From repository root, run (which will automatically download and use pre-defined conda environments):
    - `snakemake -c 1 -j 250 --use-conda --slurm --default-resources slurm_account=\<your SLURM account\> slurm_partition=\<your SLURM partition\>`, or
    - `snakemake -c 1 -j 250 --use-conda --cluster "sbatch -c {resources.cpus_per_task} --mem={resources.mem_mb}MB --time={resources.runtime} --account=\<your SLURM account\> --partition=\<your SLURM partition\>"`

### Bare Metal
- Ensure all dependencies are installed on the system. For a list reference 'workflow/envs'.
- Clone github repository.
- Create configuration file per config/README.md specifications.
- From repository root, run:
    - `snakemake -c 1 -j 250`

### My Bare Metal with SLURM commands
- `conda activate sda`
- Then I run:
    - `snakemake -c 1 -j 250 -k --slurm --default-resources slurm_account=mchaisso_100 slurm_partition=qcb`, or
    - `snakemake -c 1 -j 250 -k --use-conda --cluster "sbatch -c {resources.cpus_per_task} --mem={resources.mem_mb}MB --time={resources.runtime} --account=mchaisso_100 --partition=qcb --output=slurm-logs/slurm-%j.out"`

## Salient Output File Specifications

results/G01\_dups\_\<isoform_grouping_type\>.bed
| Column | Description |
| --- | ------ |
| #chr | Gene copy's position in assembly. |
| start | ^ |
| end | ^ |
| gene | Gene name. |
| orig\_chr | Position of gene copy's original copy in assembly. |
| orig\_start | ^ |
| orig\_end | ^ |
| strand | Strand on which gene copy is on: 0 for 'Original', '1' for reverse. |
| haplotype | The haplotype of the hit. This requires 'haplotype1' or 'haplotype2' to be explicitly stated in the chromosome/scaffold/contig name. |
| p\_identity | Similarity to 'Original' gene calculated as: #matches/(#matches+#mismatches+#insertion\_events+#deletion\_events) |
| p\_accuracy | Similarity to 'Original' gene calculated as : #matches/(#matches+#mismatches+#insertions+#deletions) |
| identity | Notes whether copy is the 'Original' or a resolved 'Copy'. |
| depth | Mean gene depth over mean assembly depth. |
| depth\_stdev | Standard deviation of gene depth over mean assembly depth calculated using 100 bp bins. |
| copy\_num | Rounded depth value. |
| depth\_by\_vcf | Copy number as determined by hmm's vcf output. |

  
results/G03\_per\_gene\_counts\_\<isoform_grouping_type\>.tsv
| Column | Description |
| --- | ------ |
| gene | Gene name. |
| depthNormalizedToAsm | Mean of depth across all gene copies divided by mean assembly depth. |
| depth | Sum of depth across all gene copies. |
| measuredCov | depth rounded to nearest whole number. |
| copyCount | Number of total gene copies - resolved plus collapsed. |
| resolvedCount | Number of resolved gene copies. |

  
results/G04\_summary\_stats\_\<isoform_grouping_type\>.tsv  
Contains tab delimited summary statistics.

## Tips

If rerunning pipeline using the same input assembly and reads, retain and copy over files A01\_assembly.fasta and A04\_assembly.bam as aligning reads to assembly is the longest step.
