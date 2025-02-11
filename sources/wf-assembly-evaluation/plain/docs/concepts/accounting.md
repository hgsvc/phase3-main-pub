# File accounting and manifest creation

A vital piece of information needed to trace problems
to their source and to "just know what has been done"
is the knowledge about all files that were processed
and generated as part of the workflow. In simple terms,
it is mandatory to keep a complete record of all input,
output/result and potential reference files for each
workflow run.

In the context of the CUBI, the process of generating
this comprehensive record is termed **"(file) accounting"**.
The automatically generated accounting information is then
combined to generate the so-called **manifest file**.

The **manifest file** contains the basic info (name, source,
size, checksums) of all of the above mentioned file types, i.e.,
input, output/result and reference files. The information is
stored as a tab-separated text file (a TSV table).

The **manifest file** will be created in the snakemake
working directory, with the name potentially modified by
a user-supplied suffix (via the `--config suffix=SUFFIX`)
option:

```
# default name
wd/manifest.tsv

# name modified with suffix
wd/manifest.SUFFIX.tsv
```

## How to generate the manifest file?

For technical reasons, the file accounting only happens
during a Snakemake dry run. After executing the workflow
twice in dry run mode, all information needed for the
accounting has been collected, and the necessary jobs
(obtaining file sizes, computing checksums etc.) are known
to Snakemake and will be processed together with the regular
data analysis tasks.

The option to execute a dry run of the workflow is the
`-n` (or `--dry-run`) parameter:

```
# execute this twice
snakemake -n [...other options...] run_all

# afterwards, start the run, the manifest
# will be created automatically
snakemake [...other options...] run_all
```

## Troubleshooting: broken file accounting

After the file accounting information has been collected,
Snakemake knows what jobs need to be processed (= which
files to generate). If, however, one of the data files
is changed (say, renamed), the file accounting will break
because the original accounting information can no longer
be generated. Consider the following (simplified) example:

```
# a data file in the workflow gives rise
# to several supplementary files needed
# for the accounting
file_A.data

# the following is just accounting info
# derived using file_A.data as input
file_A.size
file_A.source
file_A.checksum
```

Assuming `file_A.data` was renamed, the following will happen:

```
# the renamed file
file_B.data

# new accounting info
file_B.size
file_B.source
file_B.checksum

# however, the following accounting info is still
# registered, but cannot be created because
# file_A.data does not exist any more
file_A.size
file_A.source
file_A.checksum
```

In this situation, Snakemake would fail with the error that there are
missing input files (here: `file_A.data` is no longer present), and files
dependent on that input thus cannot be produced. This can be fixed by
resetting and recreating all accounting files:

```
# dry run one including the config option to
# reset the accounting information
snakemake [...other options...] -n --config resetacc=True

# dry run two after the reset as usual
snakemake [...other options...] -n
```

## Troubleshooting: it just does not work

If the manifest creation just seems impossible, (i) please get in
touch; then, instead of targeting the rule `run_all` when
executing Snakemake, (ii) please target the rule `run_all_no_manifest`.
This rule excludes the manifest creation and no problems related
to that should stop you from proceeding with your analysis.
