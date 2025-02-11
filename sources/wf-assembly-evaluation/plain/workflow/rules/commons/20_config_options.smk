import collections
import dataclasses


CFG_PARAM = collections.namedtuple("CFG_PARAM", "name default help")


@dataclasses.dataclass(frozen=True, eq=False)
class ConfigOptions:
    acclock: CFG_PARAM = CFG_PARAM(
        "wait_acc_lock_secs",
        60,
        (
            "Wait X seconds before raising/failing"
            " when attempting to acquire the file lock"
            " to update the accounting files."
            " Default: 60 [expert-only parameter]"
        ),
    )
    cpu_low: CFG_PARAM = CFG_PARAM(
        "cpu_low",
        2,
        (
            "Number of CPU cores (threads) to use"
            " for small jobs. Default: 2"
        )
    )
    cpu_medium: CFG_PARAM = CFG_PARAM(
        "cpu_medium",
        4,
        (
            "Number of CPU cores (threads) to use"
            " for medium-sized jobs. Default: 4"
        )
    )
    cpu_high: CFG_PARAM = CFG_PARAM(
        "cpu_high",
        6,
        (
            "Number of CPU cores (threads) to use"
            " for large jobs. Default: 6"
        )
    )
    cpu_max: CFG_PARAM = CFG_PARAM(
        "cpu_max",
        8,
        (
            "Number of CPU cores (threads) to use"
            " for huge jobs. This number should"
            " typically be equivalent to the maximal"
            " number of CPU core available on a single"
            " server in the targeted compute"
            " infrastructure. Default: 8"
        )
    )
    devmode: CFG_PARAM = CFG_PARAM(
        "devmode",
        False,
        (
            "Set '--config devmode=True' to ignore"
            " missing default directories underneath"
            " Snakemake's working directory."
            " Default: False"
        ),
    )
    env_singularity: CFG_PARAM = CFG_PARAM(
        "env_module_singularity",
        "Singularity",
        (
            "If the Singularity executable has to"
            " be loaded via an ENV module (common"
            " on HPC infrastructure), specify the"
            " name of the module to load."
            " Default: Singularity"
        ),
    )
    omit: CFG_PARAM = CFG_PARAM(
        "omit_options",
        ["wait_acc_lock_secs", "acclock"],
        (
            "Use this listing to suppress"
            " printing the help info for certain"
            " config options that are irrelevant"
            " to end users. [Dev-only parameter]"
        )
    )
    resetacc: CFG_PARAM = CFG_PARAM(
        "resetacc",
        False,
        (
            "Set '--config resetacc=True' to reset"
            " all accounting files (= empty them)."
            " This can be used, e.g., if a file target"
            " has been renamed and, thus, the original"
            " file can never be accounted for in the"
            " pipeline manifest output. Resetting the"
            " internal accounting files clears all"
            " entries, triggering a complete rebuild"
            " of the accounting info when executing"
            " the pipeline in dry run mode twice"
            " afterwards. Default: False"
        ),
    )
    samples: CFG_PARAM = CFG_PARAM(
        "samples",
        "",
        (
            "Specify a sample sheet at the command line"
            " via '--config samples=PATH-TO-SAMPLE-SHEET.tsv'."
            " Note that the sample sheet must be a tab-separated"
            " table (text file) with file extension '.tsv'."
            " Default: <empty>"
        ),
    )
    suffix: CFG_PARAM = CFG_PARAM(
        "suffix",
        "",
        (
            "Specify a suffix to append to the pipeline"
            " manifest output file and the run config"
            " dump. The name of the sample sheet (if"
            " provided) can be used by setting the option"
            " '--config suffix=derive'. Note that the"
            " suffix will always be reduced to letters,"
            " digits and 'minus' (hyphen) after replacing"
            " all dots and underscores with minus."
            " [Example]: sample.sheet_123.tsv => sample-sheet-123"
            " [Example]: Batch-A.098 => Batch-A-098"
            " Default: <empty>"
        ),
    )
    refcon: CFG_PARAM = CFG_PARAM(
        "use_reference_container",
        False,
        (
            "Set to true to load reference data files"
            " from reference containers. Note that setting"
            " this option to true requires also setting (i)"
            " the option 'reference_container_store' to a"
            " valid path containing all reference container"
            " files; and (ii) the option "
            " 'reference_container_names' listing all"
            " containers (by name) to be used. Default: False"
        ),
    )
    refnames: CFG_PARAM = CFG_PARAM(
        "reference_container_names",
        None,
        (
            "List the names of the reference containers"
            " that provide the required reference files."
            " Default: None"
        ),
    )
    refstore: CFG_PARAM = CFG_PARAM(
        "reference_container_store",
        None,
        (
            "Specify the path on the file system to the"
            " folder containing the reference container"
            " images (*.sif files). Default: None"
        ),
    )


OPTIONS = ConfigOptions()
