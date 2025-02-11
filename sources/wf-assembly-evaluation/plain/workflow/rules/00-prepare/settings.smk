
PATH_ID_LENGTH = 8

DATA_ROOT = config.get("data_root", "/")
DATA_ROOT = pathlib.Path(DATA_ROOT).resolve(strict=True)

CONTAINER_STORE = config.get("container_store", "/")
CONTAINER_STORE = pathlib.Path(CONTAINER_STORE).resolve(strict=True)


WILDCARDS_REF_GENOMES = list(config.get("refgenomes", dict()).keys())
if not WILDCARDS_REF_GENOMES:
    logerr("No reference genomes specified in config - entry 'refgenomes' is missing?")
    raise ValueError("No reference genomes specified")

WILDCARDS_GENE_MODELS = list(config.get("gene_models", dict()).keys())

# For analyses that require a (mostly) complete reference genome,
# use this one (e.g., gene completeness, assembly of complex regions such as rDNA)
COMPLETE_REF_GENOME = config.get("complete_reference", WILDCARDS_REF_GENOMES[0])
assert COMPLETE_REF_GENOME in WILDCARDS_REF_GENOMES

# STATISTICS PARAMETERS

SEQUENCE_LENGTH_THRESHOLDS_ASSEMBLY = {
    "default": [10000, 100000, 1000000, 10000000]
}
for asm_unit, thresholds in config.get("sequence_length_thresholds_assembly", dict()).items():
    if asm_unit.startswith("setting"):
        continue
    elif asm_unit.startswith("contaminants"):
        lookup_name = asm_unit
    elif asm_unit.startswith("asm-"):
        lookup_name = asm_unit
    elif asm_unit.startswith("asm_"):
        lookup_name = asm_unit.replace("_", "-")
    else:
        lookup_name = f"asm-{asm_unit}"
    SEQUENCE_LENGTH_THRESHOLDS_ASSEMBLY[lookup_name] = thresholds
ASSEMBLY_UNITS_PLUS_CONTAM = sorted(k for k in SEQUENCE_LENGTH_THRESHOLDS_ASSEMBLY.keys() if k != "default")
ASSEMBLY_UNITS_NO_CONTAM = sorted(k for k in SEQUENCE_LENGTH_THRESHOLDS_ASSEMBLY.keys() if k not in ["default", "contaminants"])

# This is needed to assess assembly completeness in terms
# of gene model completeness. The module 50-postprocess::asm_karyo_est
# assigns a sex to these assembly components using just the contig-to-reference
# alignments plus the info which chromosomes are indicating the respective
# sex in the karyotype
ASSEMBLY_UNITS_SEX_SPECIFIC = []
for asm_unit in config.get("sex_specific_assembly_units", []):
    if asm_unit.startswith("asm_"):
        lookup_name = asm_unit.replace("_", "-")
        ASSEMBLY_UNITS_SEX_SPECIFIC.append(lookup_name)
    elif asm_unit.startswith("asm-"):
        ASSEMBLY_UNITS_SEX_SPECIFIC.append(asm_unit)
    else:
        lookup_name = f"asm-{asm_unit}"
        if lookup_name not in ASSEMBLY_UNITS_PLUS_CONTAM:
            logerr(f"Name of sex-specific assembly unit seems to be malformed: {asm_unit} / {lookup_name}")
            raise ValueError("Name of assembly units must be 'asm_NAME' (NAME all lowercase characters).")
        ASSEMBLY_UNITS_SEX_SPECIFIC.append(lookup_name)

FORCE_ANNOTATED_SAMPLE_SEX = config.get("force_annotated_sample_sex", False)
assert isinstance(FORCE_ANNOTATED_SAMPLE_SEX, bool)

# TODO: refactor this and the above into a function
ASSEMBLY_UNITS_MAIN = []
for asm_unit in config.get("main_assembly_units", []):
    if asm_unit.startswith("asm_"):
        lookup_name = asm_unit.replace("_", "-")
        ASSEMBLY_UNITS_MAIN.append(lookup_name)
    elif asm_unit.startswith("asm-"):
        ASSEMBLY_UNITS_MAIN.append(asm_unit)
    else:
        lookup_name = f"asm-{asm_unit}"
        if lookup_name not in ASSEMBLY_UNITS_PLUS_CONTAM:
            logerr(f"Name of main assembly unit seems to be malformed: {asm_unit} / {lookup_name}")
            raise ValueError("Name of assembly units must be 'asm_NAME' (NAME all lowercase characters).")
        ASSEMBLY_UNITS_MAIN.append(lookup_name)


# TOOL PARAMETERS

## mosdepth

### mosdepth global on/off switch
RUN_MOSDEPTH = config.get("run_mosdepth", False)
MOSDEPTH_PARAMS = config.get("mosdepth", dict())

### window size for assembly-to-reference
### alignments; to evaluate assembly contig
### coverage in reference genome
RUN_MOSDEPTH_ASSM_REF_COV = RUN_MOSDEPTH and config.get("run_mosdepth_assm_ref_cov", False)
MOSDEPTH_ASSM_REF_COV_WINDOW_SIZE = int(MOSDEPTH_PARAMS.get("assm_ref_cov_window_size", 10000))

MOSDEPTH_ASSM_REF_COV_MAPQ_THRESHOLDS = [f"{t:02}" for t in MOSDEPTH_PARAMS.get("assm_ref_cov_mapq_threshold", [0, 60])]
assert isinstance(MOSDEPTH_ASSM_REF_COV_MAPQ_THRESHOLDS, list)

### similar to contig-to-reference coverage,
### here read to assembly coverage
RUN_MOSDEPTH_ASSM_READ_COV = RUN_MOSDEPTH and config.get("run_mosdepth_assm_read_cov", False)
MOSDEPTH_ASSM_READ_COV_WINDOW_SIZE = int(MOSDEPTH_PARAMS.get("assm_read_cov_window_size", 1000))

MOSDEPTH_ASSM_READ_COV_MAPQ_THRESHOLDS = [f"{t:02}" for t in MOSDEPTH_PARAMS.get("assm_read_cov_mapq_threshold", [0, 60])]
assert isinstance(MOSDEPTH_ASSM_REF_COV_MAPQ_THRESHOLDS, list)

### parameters to summarize read coverage in windows

USE_READ_TYPES_WINDOW_COVERAGE = config.get("use_read_types_window_coverage", [])
GROUPS_WINDOW_READ_COVERAGE = config.get("groups_window_read_coverage", [])

if RUN_MOSDEPTH_ASSM_READ_COV:
    assert len(USE_READ_TYPES_WINDOW_COVERAGE) > 0
    assert len(GROUPS_WINDOW_READ_COVERAGE) > 0

## NCBI FCS toolkit

RUN_NCBI_FCS_ADAPTOR = config.get("run_ncbi_fcs_adaptor", False)
RUN_NCBI_FCS_GX = config.get("run_ncbi_fcs_gx", False)

NCBI_FCS_PARAMS = config.get("ncbi_fcs", dict())
NCBI_FCS_ROOT_PATH = pathlib.Path(NCBI_FCS_PARAMS.get("root_path", "/"))
NCBI_FCS_ROOT_PATH = NCBI_FCS_ROOT_PATH.resolve(strict=True)

NCBI_FCS_ADAPTOR_SIF = pathlib.Path("CONFIG.NOT-SET.NCBI_FCS_ADAPTOR_SIF")
NCBI_FCS_ADAPTOR_SCRIPT = pathlib.Path("CONFIG.NOT-SET.NCBI_FCS_ADAPTOR_SCRIPT")
NCBI_FCS_ADAPTOR_TAXONOMY = ""
if RUN_NCBI_FCS_ADAPTOR:
    NCBI_FCS_ADAPTOR_SIF = NCBI_FCS_ROOT_PATH.joinpath(
        NCBI_FCS_PARAMS["adaptor_sif"]
    )
    assert NCBI_FCS_ADAPTOR_SIF.is_file()
    NCBI_FCS_ADAPTOR_SCRIPT = NCBI_FCS_ROOT_PATH.joinpath(
        NCBI_FCS_PARAMS["adaptor_script"]
    )
    assert NCBI_FCS_ADAPTOR_SCRIPT.is_file()
    NCBI_FCS_ADAPTOR_TAXONOMY = NCBI_FCS_PARAMS["adaptor_taxonomy"]


NCBI_FCS_GX_DB_PATH = pathlib.Path("CONFIG.NOT-SET.NCBI_FCS_GX_DB_PATH")
NCBI_FCS_GX_DB_NAME = NCBI_FCS_PARAMS.get("gx_db_name", "CONFIG.NOT-SET.NCBI_FCS_GX_DB_NAME")
NCBI_FCS_GX_SIF = pathlib.Path("CONFIG.NOT-SET.NCBI_FCS_GX_SIF")
NCBI_FCS_GX_SCRIPT = pathlib.Path("CONFIG.NOT-SET.NCBI_FCS_GX_SCRIPT")
NCBI_FCS_GX_TAX_ID = 0
if RUN_NCBI_FCS_GX:
    NCBI_FCS_GX_DB_PATH = NCBI_FCS_ROOT_PATH.joinpath(
        NCBI_FCS_PARAMS["gx_db_path"]
    )
    assert NCBI_FCS_GX_DB_PATH.is_dir()
    assert not NCBI_FCS_GX_DB_NAME.startswith("CONFIG.NOT-SET")
    NCBI_FCS_GX_SIF = NCBI_FCS_ROOT_PATH.joinpath(
        NCBI_FCS_PARAMS["gx_sif"]
    )
    assert NCBI_FCS_GX_SIF.is_file()
    NCBI_FCS_GX_SCRIPT = NCBI_FCS_ROOT_PATH.joinpath(
        NCBI_FCS_PARAMS["gx_script"]
    )
    assert NCBI_FCS_GX_SCRIPT.is_file()
    NCBI_FCS_GX_TAX_ID = NCBI_FCS_PARAMS["gx_tax_id"]


### parameters for HMMER motif search

RUN_HMMER = config.get("run_hmmer", True)
HMMER_MOTIF_SEARCH = config.get("hmmer_motif_search", dict())
HMMER_MOTIF_NAMES = sorted(HMMER_MOTIF_SEARCH.keys())

### parameters for RepeatMasker annotation

RUN_REPEATMASKER = config.get("run_repeatmasker", True)
REPEATMASKER_OFFLINE_SETUP = config.get("repeatmasker_offline_setup", False)

REPEATMASKER_SPECIES = config.get("repeatmasker_species", None)
if REPEATMASKER_SPECIES is None and RUN_REPEATMASKER:
    raise ValueError("Workflow configured to run repeatmasker, but no 'species' set")


