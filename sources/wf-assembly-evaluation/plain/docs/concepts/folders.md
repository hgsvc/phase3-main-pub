# Workflow folder structure

The default layout for the Snakemake working directory, i.e.
the directory specified with the option `--directory` / `-d`
when invoking `snakemake`, is standardized to ease navigating
through the subfolders. In the following, the abbreviation
WD is used to refer to that directory.

If you followed the CUBI recommendation and created your
Snakemake WD with the `init.py` script from the workflow
repository, the WD will have been created one level above
the workflow repository location (where the `init.py`
script resides).

Inside the WD, you find several subfolders ...

## Users: continue reading here

The only subfolder relevant to you is called
`results/`. This folder contains all important
result files of the analysis, plus two
special files (`manifest*` and `run_config*`)
that you should never delete, but can ignore.
The entire content of the `results/` folder should
be saved in a secure location (with backup). Deleting
any other folder inside the WD does not affect your results,
but may cause delays if the analysis run is changed
or continued later on - so please be sure that the
analysis is complete before deleting anything.

## Developers: continue reading here

By convention, a pipeline must sort all data files
into one of the following folders:

- `proc/`: processing; intermediate/temp data files, not vital
- `results/`: everything relevant to the client/user
  - automatically contains a copy of the sample sheet,
    the manifest file, and a dump of the run config
- `global_ref/`: global reference files either loaded from
    reference containers or *ex nihilo*, i.e. references
    manually placed in the pipeline context
- `local_ref/`: reference data local to the pipeline context,
    e.g. postprocessed files generated as part of the run
- `log/`: target for Snakemake `log:` files
- `rsrc/`: target for Snakemake `benchmark:` files

There are a few other standard directories defined that you must
use if needed. See the Snakemake module `workflow/rules/common/10_constants.smk`
for a complete list.
