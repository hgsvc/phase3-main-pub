# For CUBI developers

The following information is only relevant for CUBI developers,
and is further contextualized via the guidelines listed in
the [CUBI knowledge base](https://github.com/core-unit-bioinformatics/knowledge-base).

## Developing a Snakemake workflow locally

A more detailed explanation of the workflow setup process
can be found in the [user documentation](./running.md) for
running a workflow.

In brief, ...

1. run `./init.py --dev-only` (requires Python3)
    - this will skip creating the workflow working directory and subfolders
2. activate the created Conda environment: `conda activate ./dev_env`
3. write your code, and add tests to `workflow/snaketests.smk`
4. run tests:
    - note that some tests may be expected to fail and may produce error messages
    - if Snakemake reports a successful pipeline run, then all tests have succeeded
      irrespective of log messages that look like errors
    - if you want to test the functions loading reference data from reference containers,
      you need to build the test container `test_v0.sif` and copy it into the
      working directory for the workflow test run. Refer to the
      [reference container repository](https://github.com/core-unit-bioinformatics/reference-container)
      for build instructions.

```bash
# Example: test w/o reference container
# Note: execute the workflow first in
# '--dryrun' mode twice to trigger
# (and test) the complete MANIFEST
# creation
snakemake --cores 1 \
    [--dryrun] \
    --config devmode=True \
    --directory wd/ \
    --snakefile workflow/snaketests.smk \
    run_tests

# Example: test w/ reference container;
# the container 'test_v0.sif' must exist
# in the working directory: 'wd/test_v0.sif'
# Note: execute the workflow first in
# '--dryrun' mode twice to trigger
# (and test) the complete MANIFEST
# creation
snakemake --cores 1 \
    [--dryrun] \
    --config devmode=True \
    --directory wd/ \
    --configfiles config/testing/params_refcon.yaml \
    --snakefile workflow/snaketests.smk \
    run_tests
```

Note that if there are problems related to creating the
MANIFEST, you can first check if all other tests pass
by targeting the rule `run_tests_no_manifest`.

## Code style and checks

Please read the knowledge base articles on
[naming conventions and code style/format](https://github.com/core-unit-bioinformatics/knowledge-base/wiki/Naming-and-style).

The following helpers should be executed before committing and
pushing code to shared repositories:

1. Python scripts
  - MUST: linting: `pylint <script.py>`
  - SHOULD: organize imports: `isort <script.py>`
  - MUST: code formatting: `black [--check] <script.py>`
2. Snakemake workflow
  - MUST: code formatting `snakefmt [--check] <snakefile>`
    - note: the Snakemake formatter has some issues and sometimes
      violates Python formatting rules (see here:
      [gh#20](https://github.com/core-unit-bioinformatics/template-snakemake/issues/20)).
      The output should be manually checked and obvious idiosyncracies
      fixed before committing the code.
  - MAY: linting: `snakemake --config devmode=True --lint`
    - note: the Snakemake linter hardly produces overly helpful
      hints to improve the code quality of the Snakefile.
      Give it a try and see for yourself.
3. R scripts: tools for linting/formatting are open issues
