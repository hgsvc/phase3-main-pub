# Executing a standard CUBI workflow

## Snakemake

Snakemake workflows that implement the default interface of the
CUBI can always be deployed in the same way. The following steps
outline the recommended process:

1. create a directory for the project (`project_dir/`)
2. clone the workflow repository into that directory, and
   checkout the version of the workflow you want to run
   (if applicable)
    ```bash
    project_dir/
        |
        --- workflow_repository/
                |
                --- init.py
    ```
3. run the `init.py` script from inside the workflow repository
    - **attention**: by default, the `init.py` script attempts to create
      a [Conda environment](https://docs.conda.io/projects/conda/en/latest/index.html)
      containing all necessary tools to execute the workflow
      (that is essentially Snakemake plus a few dependencies).
      This part of the setup requires a working Conda installation
      including a proper configuration of Conda to make use of the
      [`bioconda` channel](https://bioconda.github.io/).
    - if you don't want to make use of Conda environments to run
      the workflow, please make sure that all software dependencies
      listed in the [`exec_env.yaml`](../../workflow/envs/exec_env.yaml)
      environment specification are available on your system.
4. the `init.py` script will create the following directories
    ```bash
    project_dir/
        |
        --- workflow_repository/
        |       |
        |       --- init.py
        |
        --- wd/  # the working directory for Snakemake
        --- exec_env/  # the Conda execution environment
    ```
5. activate the Conda execution environment: `conda activate ./exec_env`
6. if applicable, prepare the Snakemake profile for your compute infrastructure
    - HHU-internal users: you can make use of a small utility maintained
      by the CUBI that automates creating Snakemake profiles for the compute
      infrastructure at HHU/UKD [snakemake-utils@HHU GitLab](https://git.hhu.de/cubi/snakemake-utils)
7. prepare the sample data information (the so-called "sample sheet") as a
   plain text tab-separated table ("*.tsv"). Please refer to the specific
   workflow documentation to learn what data the respective workflow needs
   as input.
8. run the workflow from inside the workflow repository as follows:
    ```bash
    snakemake -n \  # start with a dry run
        -d ../wd/ \  # the working directory as created above
        --configfiles [...] \  # parameters for the workflow
        --config samples=PATH/TO/SAMPLE-SHEET.tsv \  # recommended: use an absolute path
        --profile PATH/TO/SNAKEMAKE-PROFILE/ \
        run_all  # create all result files specified in the workflow
    ```
    - note that executing the workflow first in dry run mode is strongly
      recommended to check if the setup process worked as expected. Moreover,
      running the workflow *twice* in dry run mode is required to create
      the [manifest output](./accounting.md) of the workflow.
9. collect the important workflow output from the `results/` folder:
    ```bash
    project_dir/
        |
        --- workflow_repository/
        |       |
        |       --- init.py
        |
        --- wd/
             |
             --- results/
        --- exec_env/
    ```
    - the `results/` folder also contains other important information
      besides the analysis output. Please refer to the description of
      the [standard folder layout](./folders.md) of the working
      directory to find more information.
