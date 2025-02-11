import dataclasses
import enum
import pathlib

# Note: this module serves as single source
# of truth regarding the constant paths
# listed below. Changing this module
# affects other code:
# - init.py script

@dataclasses.dataclass(frozen=True, eq=False)
class ConstDirectories:
    proc: pathlib.Path = pathlib.Path("proc")
    results: pathlib.Path = pathlib.Path("results")
    log: pathlib.Path = pathlib.Path("log")
    rsrc: pathlib.Path = pathlib.Path("rsrc")
    cluster_log_out: pathlib.Path = pathlib.Path("log", "cluster_jobs", "out")
    cluster_log_err: pathlib.Path = pathlib.Path("log", "cluster_jobs", "err")
    global_ref: pathlib.Path = pathlib.Path("global_ref")
    local_ref: pathlib.Path = pathlib.Path("local_ref")
    scripts: pathlib.Path = pathlib.Path("scripts")
    envs: pathlib.Path = pathlib.Path("envs")
    cache: pathlib.Path = pathlib.Path("proc", ".cache")
    accounting: pathlib.Path = pathlib.Path("proc", ".accounting")
    cache_refcon: pathlib.Path = pathlib.Path("proc", ".cache", "refcon")
    _no_init: tuple = ("scripts", "envs", "cache", "accounting", "cache_refcon")


CONST_DIRS = ConstDirectories()


@dataclasses.dataclass(frozen=True, eq=False)
class ConstFiles:
    refcon_cache: pathlib.Path = pathlib.Path(
        CONST_DIRS.cache_refcon, "refcon_manifests.cache"
    )
    account_inputs: pathlib.Path = pathlib.Path(CONST_DIRS.accounting, "inputs.listing")
    account_references: pathlib.Path = pathlib.Path(
        CONST_DIRS.accounting, "references.listing"
    )
    account_results: pathlib.Path = pathlib.Path(
        CONST_DIRS.accounting, "results.listing"
    )
    # Note that both manifest and config dump may be modified via the
    # "suffix" option and the path may thus differ at runtime.
    # Yes, yes, quite idiosyncratic, but nothing is really constant in
    # Python anyway ...
    manifest: pathlib.Path = pathlib.Path(CONST_DIRS.results, "manifest.tsv")
    config_dump: pathlib.Path = pathlib.Path(CONST_DIRS.results, "run_config.yaml")


CONST_FILES = ConstFiles()


class TimeUnit(enum.Enum):
    HOUR = 1
    hour = 1
    hours = 1
    hrs = 1
    h = 1
    MINUTE = 2
    minute = 2
    minutes = 2
    min = 2
    m = 2
    SECOND = 3
    second = 3
    seconds = 3
    sec = 3
    s = 3


class MemoryUnit(enum.Enum):
    BYTE = 0
    byte = 0
    bytes = 0
    b = 0
    B = 0
    KiB = 1
    kib = 1
    kb = 1
    KB = 1
    k = 1
    K = 1
    kibibyte = 1
    MiB = 2
    mib = 2
    mb = 2
    MB = 2
    m = 2
    M = 2
    mebibyte = 2
    GiB = 3
    gib = 3
    gb = 3
    GB = 3
    g = 3
    G = 3
    gibibyte = 3
    TiB = 4
    tib = 4
    tb = 4
    TB = 4
    t = 4
    T = 4
    tebibyte = 4
