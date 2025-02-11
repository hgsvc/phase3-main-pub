# Variant-based QV value estimation

This contains the pipeline used to run additional evaluations of consensus haplotypes for the HGSVC3 paper.
See the steps below to replicate the analysis.

## How to replicate the results for PanGenie-SHAPEIT

### Step 1: Downloading input files

Download PanGenie-based consensus haplotypes:
```bat
wget -p inputs/ http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Genotyping_1kGP/20240719_shapeit-phasing/shapeit-pangenie_chm13_all_consensus-haplotypes.agc
```

Download CHM13 reference genome:
```bat
wget -P inputs/ https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
```

Download these HGSVC3 verkko assemblies and put them into ``inputs/``:
```bat
HG00096.vrk-ps-sseq.asm-hap1.fasta.gz
HG00096.vrk-ps-sseq.asm-hap2.fasta.gz
HG01114.vrk-ps-sseq.asm-hap1.fasta.gz
HG01114.vrk-ps-sseq.asm-hap2.fasta.gz
HG01596.vrk-ps-sseq.asm-hap1.fasta.gz
HG01596.vrk-ps-sseq.asm-hap2.fasta.gz
```

Download these HPRC2 assemblies (s3://human-pangenomics/submissions/6807247E-4F71-45D8-AECE-9E5813BA1D9F--verkko-v2.2.1-release2_asms/) and put them into ``inputs/``:
```bat
HG01255.assembly.haplotype1.fasta.gz
HG01255.assembly.haplotype2.fasta.gz
HG04157.assembly.haplotype1.fasta.gz
HG04157.assembly.haplotype2.fasta.gz
```
### Step 2: run snakemake pipeline
Dependencies:
* snakemake
* singularity
* conda

```bat
snakemake -j <nr_cores> --use-conda --configfile config/config-pg.yaml
```


## How to replicate the results for 1kGP consenus haplotypes

### Step 1: Downloading input files

Download HG002 Q100 assemblies (https://github.com/marbl/HG002) and put them into ``inputs/``:
```bat
hg002v1.1.mat.fasta.gz
hg002v1.1.pat.fasta.gz
```

Produce 1kGP based consensus haplotypes by running pipeline https://github.com/eblerjana/hgsvc3/tree/main/experiments/phasing. This will produce the following file: 

```bat
{results}/haplotypes/NYGC-GRCh38/all/NYGC-GRCh38_qv_consensus-haplotypes.agc
```

Download CHM13 reference genome:
```bat
wget -P inputs/ https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
```

Add HG002 and CHM13 to archive:
```bat
agc append {results}/haplotypes/NYGC-GRCh38/all/NYGC-GRCh38_qv_consensus-haplotypes.agc inputs/chm13v2.0.fa inputs/hg002v1.1.pat.fasta.gz inputs/hg002v1.1.mat.fasta.gz > inputs/NYGC-GRCh38_qv_consensus-haplotypes_chm13_HG002.agc
```


Download GRCh38 reference genome:
```bat
wget -P inputs/ http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz
```

Download these HGSVC3 verkko assemblies and put them into ``inputs/``:
```bat
HG00096.vrk-ps-sseq.asm-hap1.fasta.gz
HG00096.vrk-ps-sseq.asm-hap2.fasta.gz
HG01114.vrk-ps-sseq.asm-hap1.fasta.gz
HG01114.vrk-ps-sseq.asm-hap2.fasta.gz
HG01596.vrk-ps-sseq.asm-hap1.fasta.gz
HG01596.vrk-ps-sseq.asm-hap2.fasta.gz
NA24385.vrk-ps-sseq.asm-hap1.fasta.gz
NA24385.vrk-ps-sseq.asm-hap2.fasta.gz
```

Download these HPRC2 assemblies (s3://human-pangenomics/submissions/6807247E-4F71-45D8-AECE-9E5813BA1D9F--verkko-v2.2.1-release2_asms/) and put them into ``inputs/``:
```bat
HG01255.assembly.haplotype1.fasta.gz
HG01255.assembly.haplotype2.fasta.gz
HG04157.assembly.haplotype1.fasta.gz
HG04157.assembly.haplotype2.fasta.gz
```

### Step 2: run snakemake pipeline
Dependencies:
* snakemake
* singularity
* conda

```bat
snakemake -j <nr_cores> --use-conda --configfile config/config-1kgp.yaml
```
