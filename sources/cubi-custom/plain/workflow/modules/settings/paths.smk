import pathlib

GPFS_BASE_HILBERT = pathlib.Path("/gpfs/project/projects/medbioinf")

TOP_ROOT_DIR = GPFS_BASE_HILBERT

WORKDIR_EVAL = TOP_ROOT_DIR.joinpath("projects/assemblies/hybrids/eval/wd")
WORKDIR_ASSEMBLY = TOP_ROOT_DIR.joinpath("projects/assemblies/hybrids/verkko/wd")

BASE_SHARE_LOCATION = TOP_ROOT_DIR.joinpath("data/00_RESTRUCTURE/shares/globus/outgoing/hgsvc")

MERQURY_RESULT_ROOT_PHASED = TOP_ROOT_DIR.joinpath(
    "data/00_RESTRUCTURE/project-centric/hgsvc/processed/20240415_merqury_results/merqury"
)

# this Merqury output folder includes the stats for the Verkko unassigned
# as a separate output run; main output is for hap1/hap2 only
MERQURY_RESULT_ROOT_ALL = TOP_ROOT_DIR.joinpath(
    "data/00_RESTRUCTURE/project-centric/hgsvc/processed/20240509_UW_merqury_flagger_complete/assembly_annotations/merqury/merqury"
)

# Inspector results for both Verkko and hifiasm in subfolders
INSPECTOR_ROOT_FOLDER = BASE_SHARE_LOCATION.joinpath(
    "ebi_upload/20230926_assembly_annotations/Inspector"
).resolve(strict=True)

# Flagger results - only in assembly coordinates
FLAGGER_ROOT_FOLDER = TOP_ROOT_DIR.joinpath(
    "data/00_RESTRUCTURE/project-centric/hgsvc/processed",
    "20240509_UW_merqury_flagger_complete/assembly_annotations/flagger"
)
FLAGGER_ROOT_FOLDER.resolve(strict=True)

SEGDUP_ROOT_FOLDER = BASE_SHARE_LOCATION.joinpath(
    "ebi_upload/20230926_assembly_annotations/uwash/segdups"
)
SEGDUP_ROOT_FOLDER.resolve(strict=True)

CENTROMERE_ANNOTATION = BASE_SHARE_LOCATION.joinpath(
    "ebi_upload/20230926_assembly_annotations/cen_coords",
    "hgsvc3_verkko_v1.4.1_hifiasm_v0.19.6_nonredundant_complete_and_accurate_active_asat_HOR_arrays_v3.list"
).resolve(strict=True)

CENREGION_ANNOTATION = BASE_SHARE_LOCATION.joinpath(
    "ebi_upload/20230926_assembly_annotations/cen_coords",
    "hgsvc3_verkko_v1.4.1_hifiasm_v0.19.6_nonredundant_complete_and_accurate_centromeric_regions_v3.list"
).resolve(strict=True)

# project repo
DIR_SNAKEFILE = pathlib.Path(workflow.basedir).resolve(strict=True)
assert DIR_SNAKEFILE.name == "workflow", DIR_SNAKEFILE

# annotations folder
DIR_ANNOTATIONS = DIR_SNAKEFILE.parent.joinpath("annotations").resolve(strict=True)


### ONE MORE DATA TABLE
ILLCNV_ANNOTATION = DIR_ANNOTATIONS.joinpath("roi", "20240725_ill_specific_CNVs.flat.bed").resolve(strict=True)
####


DIR_SCRIPTS = DIR_SNAKEFILE.joinpath("scripts").resolve(strict=True)

DIR_ENVS = DIR_SNAKEFILE.joinpath("envs").resolve(strict=True)

DIR_PROC = pathlib.Path("proc")
DIR_RES = pathlib.Path("results")
DIR_RSRC = pathlib.Path("rsrc")
DIR_LOG = pathlib.Path("log")
DIR_LOCAL_REF = pathlib.Path("local_ref")
DIR_GLOBAL_REF = pathlib.Path("global_ref")

DIR_WORKING = pathlib.Path(workflow.workdir_init).resolve(strict=True)
WORKDIR = DIR_WORKING
WD = DIR_WORKING

# prep for cluster jobs
WORKDIR.joinpath("log", "cluster_jobs", "err").mkdir(exist_ok=True, parents=True)
WORKDIR.joinpath("log", "cluster_jobs", "out").mkdir(exist_ok=True, parents=True)
