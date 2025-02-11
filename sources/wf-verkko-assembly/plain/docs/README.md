# Documentation for Snakemake workflow "genome hybrid assembly"

This workflow creates state-of-the-art genome hybrid assemblies
for diploid vertebrate species. This version of the workflow was
developed for the following scenario:

- input species: human
  - successfully tested as well: muntjac
- assembler: Verkko v1.4.1
  - hifiasm v0.19.x is not yet fully integrated
- inputs:
  - long accurate reads: PacBio HiFi (Sequel-II/Revio)
    - required coverage: at least ~40X, ideally ~60X
  - long connecting reads: Oxford Nanopore ultralong (R9)
    - required coverage: ~30X ultralong (>100 kbp) reads
  - optional input for phasing:
    - trio: kmer databases created with meryl
    - HiC: HiC short reads
    - graph/node coloring for Verkko's Rukki (GFA file)
- outputs:
  - main: whole-genome assembly, potentially phased
  - main: basic (length) statistics about assembly and long reads
  - optional: a coordinate map between the homopolymer-compressed assembly graph and the linearized plain FASTA files

If you are familiar with the templated Snakemake workflow interface
of the CUBI, please proceed directly to the
[specific instructions for configuring this workflow](workflow/README.md).

Otherwise, please also read the short user documentation in the
next section to learn how to setup and run the workflow and what
`file accounting` is. The process of `file accounting` is used
to generate a `manifest file` for the workflow run that tracks
all input, output and reference files, if applicable.

## User documentation for workflow template

All standard workflows of the CUBI implement the same user
interface (or at least aim for a highly similar interface).
Hence, before [executing the workflow](concepts/running.md),
we strongly recommend reading through the documentation
that explains how we help you to keep track of your analysis
results; we refer to this concept as
[**"file accounting"**](concepts/accounting.md). This feature
of standard CUBI workflows enables the pipeline to auto-
matically create a so-called [**"manifest"** file](concepts/accounting.md)
for your analysis run.

In case of questions, please open a GitHub issue in the repository
of the workflow you are trying to execute.

## Developer documentation

Besides reading the user documentation, CUBI developers find more
information regarding standadized workflow development in the
[developer notes](concepts/developing.md). Please keep in mind
to always cross-link that information with the guidelines
published in the
[CUBI knowledge base](https://github.com/core-unit-bioinformatics/knowledge-base/wiki/).

Please raise any issues with these guidelines "close to the code",
i.e., either open an issue in the
[knowledge base repo](https://github.com/core-unit-bioinformatics/knowledge-base)
or in the affected repo for more specific cases.
