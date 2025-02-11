# Assembly Evaluation Pipeline #

This is designed as a walkthrough of the assembly_eval pipeline to generate NucFreq plots of users chosen contigs


## Setup ##

### Config file ###

The pipeline expects to have a file called:

```
asm_eval.yaml 
```

You can see asm_eval.yaml for an example. 
Or use assembly_eval_yaml_build.py to build the yaml for you.


### Input Files ###

This pipeline requires the following input files:
	* HiFi Reads and/or ONT reads 
	* assembled genome (HiCanu or HiFiAsm)
	* bed file containing regions of interest (if producing nucfreq output)

## Running Snakefile ##
This pipeline relies on the alias snakesub which can be indicated below:

```bash
alias snakesub='mkdir -p log; snakemake -p --ri --jobname "{rulename}.{jobid}" --drmaa " -V -cwd -j y -o ./log -e ./log -l h_rt={resources.hrs}:00:00 -l mfree={resources.mem}G -pe serial {threads} -w n -S /bin/bash" -w 60'
```

To run the pipeline on the cluster and generate the alignments along with the potential misassembly/collapse regions, you can execute the following:

```bash
module load miniconda/4.9.2
snakesub -s /net/eichler/vol26/7200/software/pipelines/assembly_eval/assembly_eval.smk -j {jobs}
```

## Output Options ##

This pipeline has multiple exit points: aligning all reads to contigs, producing a bed file of potential collapses and misassemblies, and nucfreq pngs (given specific regions to produce them for). 

Each of these steps has a target rule to indicate that this is the intended end point of the pipeline. 

To generate **just the alignments**, you can use the following command:

```bash
module load miniconda/4.9.2
snakesub -s /net/eichler/vol26/7200/software/pipelines/assembly_eval/assembly_eval.smk -j {jobs} all_align
```

To generate **nucFreq pngs** for given regions, you can used the following command:
```bash
module load miniconda/4.9.2
snakesub -s /net/eichler/vol26/7200/software/pipelines/assembly_eval/assembly_eval.smk -j {jobs} all_nucfreq
```

**NOTE** The user must supply a bed file with contigs and locations to generate these images. The pipeline does not have the capacity to generate these genome-wide.


## Understanding the Output ##

This pipeline generates nucFreq plots of desired contigs. A nucFreq plot displays the coverage across a given region (y axis: coverage, x axis: position) . Black dots demonstrate the coverage of the most common basepair at that position, Red demonstrate the coverage of the second most common base pair at that position. These plots are generated through the alignments of either HiFi or ONT reads to a particular samples assembly. This pipeline allows the user to input multiple samples, for each sample you can input HiFi Reads and/or ONT reads and choose between minimap2 and winnowmap for the aligner. You will also have the flexibility to include any of the familiar nucFreq options when generating the plots. 
