"""
Module with generic global settings for the workflow.
All variables defined here should (must) be
used throughout the workflow to lower the risk
of mixing up path names, creating clashes due
to typos etc.

Operations inside this module must be limited
to trivial (one-liner) statements using
only Python's standard library functions.
In other words, just by looking up the
definition line of a certain variable in this
module, it must be obvious what the variable
should be used for (what it refers to).
"""

import pathlib
import re
import sys


##################################
# Settings for compute environment
CPU_LOW = config.get(
    OPTIONS.cpu_low.name, OPTIONS.cpu_low.default
)
assert isinstance(CPU_LOW, int)

CPU_MEDIUM = config.get(
    OPTIONS.cpu_medium.name, OPTIONS.cpu_medium.default
)
CPU_MED = CPU_MEDIUM
assert isinstance(CPU_MEDIUM, int)

CPU_HIGH = config.get(
    OPTIONS.cpu_high.name, OPTIONS.cpu_high.default
)
assert isinstance(CPU_HIGH, int)

CPU_MAX = config.get(
    OPTIONS.cpu_max.name, OPTIONS.cpu_max.default
)
assert isinstance(CPU_MAX, int)

ENV_MODULE_SINGULARITY = config.get(
    OPTIONS.env_singularity.name, OPTIONS.env_singularity.default
)

##############################################
# Settings taken from Snakemake command line
# (parameters specified in a profile, not in a config)

# A run with verbose will trigger
# a few diagnostic messages to stderr
VERBOSE = workflow.verbose
assert isinstance(VERBOSE, bool)

USE_CONDA = workflow.use_conda
assert isinstance(USE_CONDA, bool)

USE_SINGULARITY = workflow.use_singularity
assert isinstance(USE_SINGULARITY, bool)

USE_ENV_MODULES = workflow.use_env_modules
assert isinstance(USE_ENV_MODULES, bool)

# special case: the --dry-run option is not accessible
# as, e.g., an attribute of the workflow object and has
# to be extracted later (see 90_staging.smk)
DRYRUN = None


#############################################
# Settings pertaining to important path/file
# information such as: workdir, name of main
# Snakefile, path to main Snakefile, path to
# repository and so on.
# For all path resolutions below, it is
# important to realize that Snakemake is
# executing inside the working directory
# (specified via "-d" option)

# Need to know this because it determines whether
# non-existing default paths will be tolerated or not
RUN_IN_DEV_MODE = config.get(
    OPTIONS.devmode.name, OPTIONS.devmode.default
)
assert isinstance(RUN_IN_DEV_MODE, bool)

### [1: paths/repo]
# start with all paths describing the repository
DIR_SNAKEFILE = pathlib.Path(workflow.basedir).resolve(strict=True)
assert DIR_SNAKEFILE.name == "workflow"  # by convention/best practices
PATH_SNAKEFILE = pathlib.Path(workflow.main_snakefile).resolve(strict=True)
NAME_SNAKEFILE = PATH_SNAKEFILE.stem
assert DIR_SNAKEFILE.samefile(PATH_SNAKEFILE.parent)

# If the name of the snakefile is not "Snakefile" and the
# developer has not set the devmode option, print a hint
# to help diagnose path resolution errors
if NAME_SNAKEFILE != "Snakefile" and not RUN_IN_DEV_MODE:
    hint_msg = "\nDEV HINT:\n"
    hint_msg += "You are probably executing a testing pipeline,"
    hint_msg += " but you did not set the config option:\n"
    hint_msg += " '--config devmode=True'\n"
    hint_msg += "This may lead to FileNotFoundErrors for the"
    hint_msg += " subfolders expected to exist"
    hint_msg += " (proc/, results/ and so on) in the working directory.\n\n"
    sys.stderr.write(hint_msg)

DIR_REPOSITORY = DIR_SNAKEFILE.parent
DIR_REPO = DIR_REPOSITORY

# use of scripts is optional, but testing script is part
# of the template (= must resolve)
DIR_SCRIPTS = DIR_SNAKEFILE.joinpath(
    CONST_DIRS.scripts
).resolve(strict=True)

DIR_ENVS = DIR_SNAKEFILE.joinpath(
    CONST_DIRS.envs
).resolve(strict=True)

### [2: paths/workdir]
# all paths describing the working directory
# plus default subfolders
DIR_WORKING = pathlib.Path(workflow.workdir_init).resolve(strict=True)
WORKDIR = DIR_WORKING

# if the workflow is executed in development mode,
# the default paths underneath the working directory
# may not exist and that is ok
WD_PATHS_MUST_RESOLVE = not RUN_IN_DEV_MODE

WD_ABSPATH_PROCESSING = CONST_DIRS.proc.resolve(strict=WD_PATHS_MUST_RESOLVE)
WD_RELPATH_PROCESSING = WD_ABSPATH_PROCESSING.relative_to(DIR_WORKING)
DIR_PROC = WD_RELPATH_PROCESSING

WD_ABSPATH_RESULTS = CONST_DIRS.results.resolve(strict=WD_PATHS_MUST_RESOLVE)
WD_RELPATH_RESULTS = WD_ABSPATH_RESULTS.relative_to(DIR_WORKING)
DIR_RES = WD_RELPATH_RESULTS

WD_ABSPATH_LOG = CONST_DIRS.log.resolve(strict=WD_PATHS_MUST_RESOLVE)
WD_RELPATH_LOG = WD_ABSPATH_LOG.relative_to(DIR_WORKING)
DIR_LOG = WD_RELPATH_LOG

WD_ABSPATH_RSRC = CONST_DIRS.rsrc.resolve(strict=WD_PATHS_MUST_RESOLVE)
WD_RELPATH_RSRC = WD_ABSPATH_RSRC.relative_to(DIR_WORKING)
DIR_RSRC = WD_RELPATH_RSRC
DIR_BENCHMARK = DIR_RSRC

WD_ABSPATH_CLUSTERLOG_OUT = CONST_DIRS.cluster_log_out.resolve(
    strict=WD_PATHS_MUST_RESOLVE
)
WD_RELPATH_CLUSTERLOG_OUT = WD_ABSPATH_CLUSTERLOG_OUT.relative_to(DIR_WORKING)
DIR_CLUSTERLOG_OUT = WD_RELPATH_CLUSTERLOG_OUT

WD_ABSPATH_CLUSTERLOG_ERR = CONST_DIRS.cluster_log_err.resolve(
    strict=WD_PATHS_MUST_RESOLVE
)
WD_RELPATH_CLUSTERLOG_ERR = WD_ABSPATH_CLUSTERLOG_ERR.relative_to(DIR_WORKING)
DIR_CLUSTERLOG_ERR = WD_RELPATH_CLUSTERLOG_ERR

WD_ABSPATH_GLOBAL_REF = CONST_DIRS.global_ref.resolve(strict=WD_PATHS_MUST_RESOLVE)
WD_RELPATH_GLOBAL_REF = WD_ABSPATH_GLOBAL_REF.relative_to(DIR_WORKING)
DIR_GLOBAL_REF = WD_RELPATH_GLOBAL_REF

WD_ABSPATH_LOCAL_REF = CONST_DIRS.local_ref.resolve(strict=WD_PATHS_MUST_RESOLVE)
WD_RELPATH_LOCAL_REF = WD_ABSPATH_LOCAL_REF.relative_to(DIR_WORKING)
DIR_LOCAL_REF = WD_RELPATH_LOCAL_REF


##########################################
# Process sample sheet (if provided) in
# conjunction with the [run] suffix option

SAMPLE_SHEET_PATH = config.get(OPTIONS.samples.name, OPTIONS.samples.default)
SAMPLE_SHEET_NAME = None

# Note, default for run suffix is empty string
RUN_SUFFIX = config.get(OPTIONS.suffix.name, OPTIONS.suffix.default)
if RUN_SUFFIX == "derive" and not SAMPLE_SHEET_PATH:
    raise ValueError("No sample sheet specified, hence, cannot derive run suffix.")

if SAMPLE_SHEET_PATH:
    SAMPLE_SHEET_PATH = pathlib.Path(SAMPLE_SHEET_PATH).resolve(strict=True)
    assert SAMPLE_SHEET_PATH.name.lower().endswith(".tsv"), \
        "Only TSV tables allowed as sample sheet (*.tsv)."
    SAMPLE_SHEET_NAME = SAMPLE_SHEET_PATH.stem
    assert SAMPLE_SHEET_NAME
    if RUN_SUFFIX == "derive":
        RUN_SUFFIX = SAMPLE_SHEET_NAME
    # set path in results/ folder to keep a copy
    # of the sample sheet as part of the output
    COPY_SAMPLE_SHEET_RELPATH = DIR_RES.joinpath(f"{SAMPLE_SHEET_NAME}.tsv")
    COPY_SAMPLE_SHEET_ABSPATH = COPY_SAMPLE_SHEET_RELPATH.resolve()
else:
    # set target path of sample sheet under results/
    # to empty path if not sample sheet provided
    COPY_SAMPLE_SHEET_RELPATH = "no-sample-sheet"
    COPY_SAMPLE_SHEET_ABSPATH = ""

# Postprocess the run suffix to consist
# only of digits, chars, and "minus"
RUN_SUFFIX = RUN_SUFFIX.replace(".", "-").replace("_", "-")
RUN_SUFFIX = "".join(re.findall("[a-z0-9\-]+", RUN_SUFFIX, re.IGNORECASE))
# in case the above resulted in two or more
# consecutive hyphens, replace with single one
RUN_SUFFIX = re.sub("\-\-+", "-", RUN_SUFFIX)

if RUN_SUFFIX:
    RUN_SUFFIX = f".{RUN_SUFFIX}"


##############################
# Settings of fixed file names
# (modulo the run suffix) for
# manifest and config dump

RUN_CONFIG_ABSPATH = DIR_WORKING.joinpath(
    CONST_FILES.config_dump.with_suffix(f"{RUN_SUFFIX}.yaml")
).resolve(strict=False)
RUN_CONFIG_RELPATH = RUN_CONFIG_ABSPATH.relative_to(DIR_WORKING)

MANIFEST_ABSPATH = DIR_WORKING.joinpath(
    CONST_FILES.manifest.with_suffix(f"{RUN_SUFFIX}.tsv")
).resolve(strict=False)
MANIFEST_RELPATH = MANIFEST_ABSPATH.relative_to(DIR_WORKING)


#############################
# Settings pertaining to the
# file accounting needed for
# creating the manifest file

# For file-based locking of accounting files.
# Expert-only option in case accounting outside
# of dry runs becomes necessary
WAIT_ACC_LOCK_SECS = config.get(
    OPTIONS.acclock.name, OPTIONS.acclock.default
)

# should the accounting files be reset/emptied?
RESET_ACCOUNTING = config.get(OPTIONS.resetacc.name, OPTIONS.resetacc.default)
assert isinstance(RESET_ACCOUNTING, bool)

# this is kept for backward-compatibility, but
# could be replaced with something more elegant
# making use of the new "10_constants" module
ACCOUNTING_FILES = {
    CONST_FILES.account_inputs.stem: CONST_FILES.account_inputs,
    CONST_FILES.account_references.stem: CONST_FILES.account_references,
    CONST_FILES.account_results.stem: CONST_FILES.account_results,
}


#####################################################
# Settings related to the use of reference containers

USE_REFERENCE_CONTAINER = config.get(
    OPTIONS.refcon.name, OPTIONS.refcon.default
)
USE_REFCON = USE_REFERENCE_CONTAINER  # shorthand
assert isinstance(USE_REFERENCE_CONTAINER, bool)

if USE_REFERENCE_CONTAINER:
    try:
        DIR_REFERENCE_CONTAINER = config[OPTIONS.refstore.name]
    except KeyError:
        raise KeyError(
            "The config option 'use_reference_container' is set to True. "
            "Consequently, the option 'reference_container_store' must be "
            "set to an existing folder on the file system containing the "
            "reference container images (*.sif files)."
        )
    else:
        DIR_REFERENCE_CONTAINER = pathlib.Path(DIR_REFERENCE_CONTAINER).resolve(
            strict=True
        )
else:
    DIR_REFERENCE_CONTAINER = pathlib.Path("/")
DIR_REFCON = DIR_REFERENCE_CONTAINER  # shorthand
