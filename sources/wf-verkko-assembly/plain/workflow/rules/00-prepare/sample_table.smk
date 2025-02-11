import pathlib
import hashlib
import collections
import re

import pandas

SAMPLES = None
MAP_SAMPLE_TO_INPUT_FILES = None

TRIO_SAMPLES = None
SSEQ_SAMPLES = None
HIC_SAMPLES = None
UNPHASED_SAMPLES = None

CONSTRAINT_TRIO_SAMPLES = None
CONSTRAINT_SSEQ_SAMPLES = None
CONSTRAINT_HIC_SAMPLES = None
CONSTRAINT_UNPHASED_SAMPLES = None


def process_sample_sheet():

    SAMPLE_SHEET_FILE = pathlib.Path(config["samples"]).resolve(strict=True)
    SAMPLE_SHEET = pandas.read_csv(
        SAMPLE_SHEET_FILE,
        sep="\t",
        header=0,
        comment="#"
    )

    # step 1: each row is a sample,
    # just collect the input files
    sample_input, sample_sets = collect_input_files(SAMPLE_SHEET)
    all_samples = sorted(sample_input.keys())
    assert len(all_samples) == SAMPLE_SHEET.shape[0]

    global SAMPLES
    SAMPLES = all_samples
    global TRIO_SAMPLES
    TRIO_SAMPLES = sample_sets["trio"]
    global SSEQ_SAMPLES
    SSEQ_SAMPLES = sample_sets["sseq"]
    global HIC_SAMPLES
    HIC_SAMPLES = sample_sets["hic"]
    global UNPHASED_SAMPLES
    UNPHASED_SAMPLES = sample_sets["unphased"]

    global MAP_SAMPLE_TO_INPUT_FILES
    MAP_SAMPLE_TO_INPUT_FILES = sample_input

    global CONSTRAINT_TRIO_SAMPLES
    CONSTRAINT_TRIO_SAMPLES = _build_constraint(TRIO_SAMPLES)
    global CONSTRAINT_SSEQ_SAMPLES
    CONSTRAINT_SSEQ_SAMPLES = _build_constraint(SSEQ_SAMPLES)
    global CONSTRAINT_HIC_SAMPLES
    CONSTRAINT_HIC_SAMPLES = _build_constraint(HIC_SAMPLES)
    global CONSTRAINT_UNPHASED_SAMPLES
    CONSTRAINT_UNPHASED_SAMPLES = _build_constraint(UNPHASED_SAMPLES)

    return


def collect_input_files(sample_sheet):
    """
    The output of this function should
    be sufficient to run the workflow
    in single-sample mode
    """
    sample_input = collections.defaultdict(dict)
    trio_samples = set()
    sseq_samples = set()
    hic_samples = set()
    unphased_samples = set()

    for row in sample_sheet.itertuples():
        sample = row.sample
        hifi_input, hifi_hashes = collect_sequence_input(row.hifi)
        ont_input, ont_hashes = collect_sequence_input(row.ont)
        if row.target == "trio":
            assert hasattr(row, "hap1") and hasattr(row, "hap2"), (
                f"Trio-phasing a sample requires hap1/hap2 "
                "columns in sample sheet: {sample}"
            )
            trio_samples.add(sample)
            hap1_db = row.hap1
            if hap1_db.endswith("meryl"):
                assert pathlib.Path(hap1_db).resolve(strict=True).is_dir()
                hap1_input_db = hap1_db
            elif hap1_db.endswith("tar.gz"):
                assert pathlib.Path(hap1_db).resolve(strict=True).is_file()
                # this path: see module in
                # 00-prepare::verkko.smk
                hap1_input_db = str(DIR_PROC.joinpath(
                    "00-prepare", "verkko", "hapmer_dbs",
                    f"{sample}.hap1.meryl"))
                sample_input[sample]["hap1_tar"] = hap1_db
            else:
                raise ValueError(f"Cannot process k-mer / hap-mer meryl DB input: {hap1_db}")
            sample_input[sample]["hap1"] = hap1_input_db

            hap2_db = row.hap2
            if hap2_db.endswith("meryl"):
                assert pathlib.Path(hap2_db).resolve(strict=True).is_dir()
                hap2_input_db = hap2_db
            elif hap2_db.endswith("tar.gz"):
                assert pathlib.Path(hap2_db).resolve(strict=True).is_file()
                # this path: see module in
                # 00-prepare::verkko.smk
                hap2_input_db = str(DIR_PROC.joinpath(
                    "00-prepare", "verkko", "hapmer_dbs",
                    f"{sample}.hap2.meryl"))
                sample_input[sample]["hap2_tar"] = hap2_db
            else:
                raise ValueError(f"Cannot process k-mer / hap-mer meryl DB input: {hap2_db}")
            sample_input[sample]["hap2"] = hap2_input_db

            assert sample_input[sample]["hap1"] != sample_input[sample]["hap2"]

        elif row.target == "sseq":
            # NB: at the moment, phasing an assembly with
            # Strand-seq requires to complete the unphased assembly
            # first, hence the sample is always added to the
            # unphased_samples set as well. This may need to be
            # changed if both pipelines can be integrated.
            unphased_samples.add(sample)
            sseq_samples.add(sample)
            assert hasattr(row, "phasing_paths"), (
                "Strand-seq phasing of sample requires"
                f"phasing_paths column in sample sheet: {sample}"
            )
            phasing_paths_file = pathlib.Path(row.phasing_paths).resolve(strict=True)
            assert phasing_paths_file.is_file(), (
                f"Phasing paths entry is not a valid file: {phasing_paths_file}"
            )
            assert phasing_paths_file.suffix == ".gaf", (
                f"Phasing paths file must be GAF: {phasing_paths_file.name}"
            )
            sample_input[sample]["phasing_paths"] = phasing_paths_file

        elif row.target == "hic":
            assert hasattr(row, "hic1"), "Phasing with Hi-C requires field 'hic1' in sample table"
            assert hasattr(row, "hic2"), "Phasing with Hi-C requires field 'hic2' in sample table"
            hic1_read_files, hic2_file_hashes = collect_sequence_input(row.hic1)
            hic2_read_files, hic2_file_hashes = collect_sequence_input(row.hic2)
            sample_input[sample]["hic1"] = hic1_read_files
            sample_input[sample]["hic2"] = hic2_read_files
            hic_samples.add(sample)

        elif row.target == "unphased":
            unphased_samples.add(sample)
        else:
            raise ValueError(row.target)
        sample_input[sample]["hifi"] = hifi_input
        sample_input[sample]["ont"] = ont_input

    sample_sets = {
        "trio": trio_samples,
        "sseq": sseq_samples,
        "hic": hic_samples,
        "unphased": unphased_samples
    }

    return sample_input, sample_sets


def _read_input_files_from_fofn(fofn_path):
    """Read input file listing from
    file of file names
    TODO: candidate for inclusion in template
    """

    input_files = []
    with open(fofn_path, "r") as listing:
        for line in listing:
            if not line.strip():
                continue
            try:
                file_path = pathlib.Path(line.strip()).resolve(strict=True)
            except FileNotFoundError:
                try:
                    file_path = DATA_ROOT.joinpath(line.strip()).resolve(strict=True)
                except FileNotFoundError:
                    err_msg = "\nERROR\n"
                    err_msg += f"Cannot find file: {line.strip()}\n"
                    err_msg += f"Data root is set to: {DATA_ROOT}\n"
                    sys.stderr.write(err_msg)
                    raise
            input_files.append(file_path)

    return sorted(input_files)


def subset_path(full_path):
    """This helper exists to reduce
    the absolute path to a file
    to just the file name and its
    parent.
    TODO: should be codified as part
    of the template utilities to improve
    infrastructure portability of active
    workflows
    """
    folder_name = full_path.parent.name
    file_name = full_path.name
    subset_path = f"{folder_name}/{file_name}"
    # if it so happens that the file resides
    # in a root-level location, strip off
    # leading slash
    return subset_path.strip("/")


def collect_sequence_input(path_spec):
    """
    Generic function to collect HiFi or ONT/Nanopore
    input (read) files
    """
    input_files = []
    input_hashes = []
    for sub_input in path_spec.split(","):
        input_path = pathlib.Path(sub_input).resolve()
        if input_path.is_file() and input_path.name.endswith(".fofn"):
            fofn_files = _read_input_files_from_fofn(input_path)
            fofn_hashes = [
                hashlib.sha256(
                    subset_path(fp).encode("utf-8")
                ).hexdigest() for fp in fofn_files
            ]
            input_files.extend(fofn_files)
            input_hashes.extend(fofn_hashes)
        elif input_path.is_file():
            input_hash = hashlib.sha256(
                subset_path(input_path).encode("utf-8")
            ).hexdigest()
            input_files.append(input_path)
            input_hashes.append(input_hash)
        elif input_path.is_dir():
            collected_files = _collect_files(input_path)
            collected_hashes = [
                hashlib.sha256(
                    subset_path(fp).encode("utf-8")
                ).hexdigest() for fp in collected_files
            ]
            input_files.extend(collected_files)
            input_hashes.extend(collected_hashes)
        else:
            raise ValueError(f"Cannot handle input: {sub_input}")
    return input_files, input_hashes


def _build_constraint(values):
    escaped_values = sorted(map(re.escape, map(str, values)))
    constraint = "(" + "|".join(escaped_values) + ")"
    return constraint


def _collect_files(folder):

    all_files = set()
    for pattern in config["input_file_ext"]:
        pattern_files = set(folder.glob(f"**/*.{pattern}"))
        all_files = all_files.union(pattern_files)
    all_files = [f for f in sorted(all_files) if f.is_file()]
    if len(all_files) < 1:
        raise ValueError(f"No input files found underneath {folder}")
    return all_files

process_sample_sheet()
