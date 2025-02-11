# CDR-Finder
This repository contains a snakemake to identify and quantify hypomethylated regions within centromeres, or Centromere Dip Regions (CDRs; Altemose et al., Science, 2022).

Original pipeline constructed by @arozanski97 with help from Glennis Logsdon

Adapted and distributed by @fkmastrorosa and @wharvey31

This is done by:
- Extracting the sequence of interest
- Intersect the sequence with its methylation data
- Bin the region of interest into 5 kbp windows and calculate their mean methylation percentage
- Run RepeatMasker on the sequence of interest to identify regions containing Alpha-satellite (ALR/Alpha)
- For each Alpha-satellite containing region, it identifies bins with a lower methylation percentage than the avearge of the region
- Merge consecutive bins
- Checks if flanking bins with greater than maximum methylation percentage (defined as within 1 standard deviation from the maximum)
- Checks if the CDR calls are smaller than a user-specified threshold

# Input
- fasta: sample genome assembly
- target_bed: BED file of target region coordinates
- meth_tsv: modbam2bed methylation BED

# Output
Based on the mean methylation frequency across the region, identifies and bins regions with methylation frequency below mean spanning >25 kbp.
```
results/{sample}_CDR.bed
```

# Requirements
This pipeline requires snakemake, singularity, and the python packages pandas and numpy.
