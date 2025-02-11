# MHC-annotation
Tools to annotate haplotypes of MHC with gene and transcript information 

## Dependencies
There is only a few dependencies which you have to take care of yourself. 

### conda environment
The easiest way to take of these dependencies is to utilize the environment.yml in this repository to create a conda environment.
`conda create --name mhca -f environment_mhca.yml`

### manually
The script assumes that you have **minimap2** installed and that it is callable with `minimap2` on the command line.
There are several ways to install this dependency depending on your environment. Giving two examples:
- `sudo apt-get install minimap2`
- `conda install -c bioconda minimap2`

You will also need the python packages *biopython* and *bcbio-gff*. You can install them in a conda environment with:

- `conda install -c bioconda bcbio-gff`
- `conda install -c bioconda biopython`

## Install
Install this package with pip.

`pip install MHC-annotation`

## Usage

This tools assumes that you are trying to annotate a sequence that stems from the human MHC region (also known as: HLA-region). If you supply any other sequence, you will most likely receive no meaningful result.
You can annotate a fasta file with the following command:

`mhca annotate <your_sequence.fa> <output_folder>`

You can check whether the transcript in the resulting gff file would produce meaningful products.

`mhca check_CDS <your_sequence.fa> <output_folder/your_sequence.gff>`
