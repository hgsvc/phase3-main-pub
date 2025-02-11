"""
Use this module to extend the default
workflow output (a list of target files)
per sub-module.
The WORKFLOW_OUTPUT list is referenced
in the main Snakefile
"""

WORKFLOW_OUTPUT = []

WORKFLOW_OUTPUT.extend(READ_STATS)

WORKFLOW_OUTPUT.extend(VERKKO_OUTPUT)
