import pathlib

def collect_merqury_kmer_tracks(wildcards):

    sample_id, assm = wildcards.sample.split(".")
    assembler = {"vrk-ps-sseq": "verkko", "hsm-ps-sseq": "hifiasm"}[assm]

    kmer_tracks = sorted(MERQURY_RESULT_ROOT_ALL.joinpath(
        f"{sample_id}-{assembler}"
        ).glob("**/*only.bed")
    )
    if assembler == "verkko":
        kmer_tracks += sorted(MERQURY_RESULT_ROOT_ALL.joinpath(
            f"{sample_id}-{assembler}-unassigned"
            ).glob("**/*only.bed")
        )

    files = []
    for kmer_track in kmer_tracks:
        assert assembler in kmer_track.name or assm in kmer_track.name, [p.name for p in kmer_tracks]
        files.append(kmer_track)
    if assembler == "verkko":
        if not len(files) == 3:
            raise RuntimeError(f"Merqury verkko error: {wildcards}")
    if assembler == "hifiasm":
        if not len(files) == 2:
            raise RuntimeError(f"Merqury hifiasm error: {wildcards}")
    return sorted(files)


def get_assessem_cli_parameters(input_files, get_list):

    track_files = []
    track_labels = []
    score_columns = []
    out_list = None
    for input_file in input_files:
        file_name = pathlib.Path(input_file).name
        if file_name.endswith(".sizes.txt"):
            continue
        track_files.append(input_file)
        if file_name.endswith(".flagger-labels.tsv.gz"):
            track_labels.append("flagger")
            score_columns.append("score")
        elif file_name.endswith(".flagger-binary.tsv.gz"):
            track_labels.append("flagbin")
            score_columns.append("binary")
        elif file_name.endswith(".hifi.inspector-errors.tsv.gz"):
            track_labels.append("inspect_hifi")
            score_columns.append("binary")
        elif file_name.endswith(".ont.inspector-errors.tsv.gz"):
            track_labels.append("inspect_ont")
            score_columns.append("binary")
        elif file_name.endswith(".merqury-asmonly-kmer.tsv.gz"):
            track_labels.append("merqury")
            score_columns.append("binary")
        elif file_name.endswith(".mosdepth-windowed.tsv.gz"):
            # NA18989.vrk-ps-sseq.hifi.mq00.mosdepth-windowed.tsv.gz
            parts = file_name.split(".")
            read_type = parts[-5]
            assert read_type in ["hifi", "ont"]
            mapq = parts[-4]
            track_labels.append(f"rd_{read_type}_{mapq}")
            score_columns.append("cov")
        elif file_name.endswith(".nucfreq-cov-bin.tsv.gz"):
            track_labels.append("nucfreq")
            score_columns.append("score")
        elif file_name.endswith(".nucfreq-binary.tsv.gz"):
            track_labels.append("nfbin")
            score_columns.append("binary")
        elif file_name.endswith(".sd-095.tsv.gz"):
            track_labels.append("sd95")
            score_columns.append("binary")
        elif file_name.endswith(".sd-098.tsv.gz"):
            track_labels.append("sd98")
            score_columns.append("binary")
        elif file_name.endswith("sseq-switch-breaks.bed"):
            track_labels.append("sseqbrkp")
            score_columns.append("score")
        elif file_name.endswith("active_asat_HOR_arrays_v3.bed"):
            track_labels.append("centro")
            score_columns.append("score")
        elif file_name.endswith("busco-issues.bed"):
            track_labels.append("busco")
            score_columns.append("score")
        else:
            raise ValueError(f"Cannot process filename: {file_name}")
    if get_list == "files":
        out_list = track_files
    if get_list == "labels":
        out_list = track_labels
    if get_list == "columns":
        out_list = score_columns
    assert out_list is not None
    return out_list


def binsize_to_int(bin_size):

    if bin_size.endswith("k"):
        factor = int(1e3)
        base_number = int(bin_size.strip("k"))
    else:
        try:
            base_number = int(bin_size)
            factor = 1
        except ValueError:
            raise ValueError(f"Cannot normalize bin size: {bin_size}")
    num_bin_size = int(base_number * factor)
    return num_bin_size


def get_assessem_mem(bin_size):

    num_bin_size = binsize_to_int(bin_size)

    if num_bin_size > int(5e4):
        mem_est = 16384
    elif num_bin_size > int(9e3):
        mem_est = 24576
    elif num_bin_size > int(5e2):
        mem_est = 32768
    else:
        mem_est = 163840
    return mem_est


def get_assessem_hrs(bin_size):

    num_bin_size = binsize_to_int(bin_size)

    if num_bin_size > int(5e4):
        hrs_est = 6
    elif num_bin_size > int(9e3):
        hrs_est = 36
    else:
        hrs_est = 48
    return hrs_est


def determine_embed_feature_set(wildcard_fset):
    """Known features for embedding:
    busco, centro, flagger, inspect_hifi, inspect_ont, merqury,
    nucfreq, rd_hifi_mq00, rd_hifi_mq60, rd_ont_mq00, rd_ont_mq60,
    sd95, sd98, sseqbrkp
    """

    defined_sets = {
        "full": "--skip-features nfbin flagbin",
        "no-mq0": "--skip-features rd_hifi_mq00 rd_ont_mq00 nfbin flagbin",
        "no-rd": "--skip-features rd_hifi_mq00 rd_hifi_mq60 rd_ont_mq00 rd_ont_mq60 nfbin flagbin",
        "mrgfull": "--skip-features nucfreq flagger --merge-features nfbin flagbin merqury inspect_hifi inspect_ont"
    }

    changed_setting = defined_sets[wildcard_fset]

    return changed_setting
