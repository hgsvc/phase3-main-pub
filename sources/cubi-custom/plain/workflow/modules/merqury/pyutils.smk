import pathlib


def find_merqury_output_file(sample, which, phased_only=True):

    sample_id = sample.split(".")[0]
    if "vrk" in sample:
        assert ASSEMBLER == "verkko"
        subfolder = f"{sample_id}-verkko"
        base_pattern = [
            f"{sample}",
            f"{sample_id}-verkko"
        ]
    elif "hsm" in sample:
        assert ASSEMBLER == "hifiasm"
        subfolder = f"{sample_id}-hifiasm"
        base_pattern = [
            f"{sample}",
            f"{sample_id}-hifiasm"
        ]
    else:
        raise NotImplementedError(sample)

    if not phased_only:
        assert ASSEMBLER == "verkko"
        search_folder = MERQURY_RESULT_ROOT_ALL.joinpath(subfolder).resolve(strict=True)
    else:
        search_folder = MERQURY_RESULT_ROOT_PHASED.joinpath(subfolder).resolve(strict=True)

    selected_file = []
    if which == "completeness":
        if phased_only:
            candidate = search_folder.joinpath("completeness.stats").resolve(strict=True)
        else:
            candidate = search_folder.joinpath(f"{sample_id}-verkko.completeness.stats").resolve(strict=True)
        selected_file.append(candidate)
    elif which == "qv_detail":
        candidates = list(search_folder.glob("*hap*.qv"))
        if ASSEMBLER == "hifiasm":
            candidates = [c for c in candidates if "-hifiasm-" in c.name]
        if ASSEMBLER == "verkko":
            candidates = [c for c in candidates if ".vrk-ps-sseq.asm-" in c.name]
        if len(candidates) == 2:
            selected_file = candidates
    elif which == "qv_summary":
        for pattern in base_pattern:
            file_candidate = search_folder.joinpath(f"{pattern}.qv")
            if not file_candidate.is_file():
                continue
            selected_file.append(file_candidate)
            break
    elif which == "errors":
        candidates = list(search_folder.glob("*_only.bed"))
        if ASSEMBLER == "verkko":
            candidates = [c for c in candidates if ".vrk-ps-sseq.asm-" in c.name]
        if ASSEMBLER == "hifiasm":
            candidates = [c for c in candidates if "-hifiasm-" in c.name]
        if len(candidates) == 2:
            selected_file = candidates
    else:
        raise NotImplementedError(which)

    if len(selected_file) < 1:
        raise FileNotFoundError(f"{sample} / {which}")
    assert isinstance(selected_file, list)
    return sorted(selected_file)


# to keep track of which input file was used
# exactly, most post-processed Merqury outputs
# will contain a file listing in the header
def shorten_merqury_file_path(file_path):

    cut_parent = pathlib.Path(file_path).parents[2]
    assert "merqury" in cut_parent.name
    short_path = str(file_path).replace(str(cut_parent), "")
    return short_path
