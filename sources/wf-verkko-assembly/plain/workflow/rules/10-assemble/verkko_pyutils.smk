
def assemble_verkko_screen_string():
    """This helper function exists to
    properly assemble the command line
    string for Verkko's >screen< option,
    which only has a default for human.
    """

    screen_files_exist = config.get("verkko_screen_files_exist", True)

    screen_opt = ""
    screen_spec = config.get("verkko_screen", None)

    if screen_spec is None:
        screen_opt = ""
    elif screen_spec == "human":
        screen_opt = " --screen human "
    elif isinstance(screen_spec, list):
        assert isinstance(screen_spec[0], dict)
        screen_opt = " "
        for spec in screen_spec:
            label, fasta = [(k, v) for k, v in spec.items()][0]
            fasta = pathlib.Path(fasta).resolve(strict=screen_files_exist)
            screen_opt += f"--screen {label} {fasta} "
    elif isinstance(screen_spec, dict):
        screen_opt = " "
        for label, fasta in screen_spec.items():
            fasta = pathlib.Path(fasta).resolve(strict=screen_files_exist)
            screen_opt += f"--screen {label} {fasta} "
    else:
        raise ValueError(f"Cannot parse Verkko screen specification: {screen_spec}")

    if VERBOSE:
        logout(f"Verkko screen option set to: {screen_opt}")

    return screen_opt


def increase_mbg_resources(attempt):
    """For some HiFi datasets, MBG requires
    much more time to build the initial graph.
    This helper function exists to directly
    increase the MBG resource requirements
    if the Verkko run is restarted.

    Note that these values are scaled again internally
    by Verkko, see Snakefiles::functions.sm::getMBGMemoryRequest

    In Verkko v1.3.1, this pertains to ...
    ... option string `--mbg-run`
    and affects settings ...
    #  build-graph
    mbg_n_cpus:          '4'
    mbg_mem_gb:          '0'
    mbg_time_h:          '72'
    ... for rule `1-buildGraph.sm::buildGraph`
    """
    mbg_resources = ""
    if int(attempt) == 2:
        # this is CPU - MEM_GB - TIME_HRS
        mbg_resources = "--mbg-run 8 160 72"
    if int(attempt) > 2:
        # this is CPU - MEM_GB - TIME_HRS
        mbg_resources = "--mbg-run 8 224 96"
    return mbg_resources


def increase_process_ont_resources(attempt):
    """Singular cases where processing the ONT
    paths is OOM killed. Use the same strategy
    as above and increase the initial job
    resources if Verkko is restarted
    (NB: Verkko also restarts its jobs)

    In Verkko v1.3.1, this pertains to...
    ... option string `--pop-run`
    and affects settings ...
    #  process_ont_paths
    pop_n_cpus:          '1'
    pop_mem_gb:          '64'
    pop_time_h:          '24'
    ... for rule `4-processONT.sm::processONT`
    """
    pop_resources = ""
    if int(attempt) > 1:
        # this is CPU - MEM_GB - TIME_HRS
        pop_resources = "--pop-run 1 96 48"
    return pop_resources


def increase_layout_contigs_resources(attempt):
    """Singular cases where creating the layout
    is OOM killed. Use the same strategy
    as above and increase the initial job
    resources if Verkko is restarted
    (NB: Verkko also restarts its jobs)

    In Verkko v1.3.1, this pertains to...
    ... option string `--lay-run`
    and affects settings ...
    #  process_ont_paths
    lay_n_cpus=1
    lay_mem_gb=32
    lay_time_h=24
    ... for rule `6-layoutContigs.sm::layoutContigs`
    """
    lay_resources = ""
    if int(attempt) > 1:
        # this is CPU - MEM_GB - TIME_HRS
        lay_resources = "--lay-run 1 96 24"
    return lay_resources


def get_verkko_output(file_collection, which_file, relpath=True):

    _this_func = "10-assemble::verkko_pyutils.smk::get_verkko_output"

    get_path = "rel_path" if relpath else "abs_path"

    if not pathlib.Path(file_collection).is_file():
        output_file_path = "file-collection-not-yet-created"
    else:
        import json

        with open(file_collection, "r") as dump:
            output_files = json.load(dump)
            try:
                output_file_path = output_files[which_file][get_path]
            except KeyError:
                err_msg = (
                    f"\nERROR in {_this_func}\n"
                    f"Cannot retrieve file '{which_file}' "
                    "from Verkko output file collection:\n"
                    f"{file_collection}\n"
                )
                logerr(err_msg)
                raise

    return output_file_path


def get_verkko_asm_units(phasing_state):
    """Borderline: no need to screen for all three
    (rdna, ebv, mito), so this will fail in those
    cases. Not obvious how to avoid that.
    TODO fix the above
    """
    if VERKKO_SCREEN:
        scraps = ["disconnected", "rdna", "ebv", "mito"]
    else:
        scraps = ["disconnected"]

    if phasing_state in ["ps-sseq", "sseq", "ps-trio", "hic", "ps-hic"]:
        asm_units = ["hap1", "hap2", "unassigned"] + scraps
    elif phasing_state in ["ps-none", "none"]:
        # Important: the whole-genome assembly FASTA represents the whole
        # genome in the Verkko space, but not in the biological space,
        # i.e. scrap sequences (rDNA, mito etc.) are not included in that file.
        asm_units = ["wg"] + scraps
    else:
        raise ValueError(f"Unknown phasing state: {phasing_state}")

    return asm_units
