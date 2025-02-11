# Constants such as standard
# folders inside the working
# directory
include: "10_constants.smk"
# Description and default
# values for all supported
# config options
include: "20_config_options.smk"
# Module defining a large number
# of global variables such as
# relative paths to be used
# throughout the workflow
include: "30_settings.smk"
# Module containing Python-only
# helper functions
include: "40_pyutils.smk"
# Module containing generic
# Snakemake rules to accomplish
# tasks such as creating the
# manifest file
include: "50_smkutils.smk"
# Module containing Snakemake
# rules concerned with dealing
# with reference containers
include: "70_refcon.smk"
# Module performing state/env-altering
# operations before Snakemake starts its
# actual work
include: "90_staging.smk"
