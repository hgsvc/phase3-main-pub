
def get_repeatmasker_run_memory_mb(input_size_mb, compressed=False):

    threshold_tiny = 1
    threshold_small = 100
    threshold_normal = 1000

    if compressed:
        compression_scaling = 1
    else:
        compression_scaling = 4  # ~gzip compressed FASTA vs uncompressed

    if input_size_mb < threshold_tiny * compression_scaling:
        mem_mb = 16384
    elif input_size_mb < threshold_small * compression_scaling:
        mem_mb = 24576
    elif input_size_mb < threshold_normal * compression_scaling:
        mem_mb = 32768
    else:
        mem_mb = 229376
    return mem_mb


def get_repeatmasker_run_time_hrs(input_size_mb, compressed=False):

    threshold_tiny = 1
    threshold_small = 100
    threshold_normal = 1000

    if compressed:
        compression_scaling = 1
    else:
        compression_scaling = 4  # ~gzip compressed FASTA vs uncompressed

    if input_size_mb < threshold_tiny * compression_scaling:
        time_hrs = 0
    elif input_size_mb < threshold_small * compression_scaling:
        time_hrs = 1
    elif input_size_mb < threshold_normal * compression_scaling:
        time_hrs = 71
    else:
        time_hrs = 84
    return time_hrs


def get_num_threads_hmmer(motif_name):
    """This separate helper function exists
    to encapsulate the THREADS value set for
    a default run, i.e. CPU_LOW
    """
    num_threads = int(min(CPU_MAX, CPU_LOW * hmmer_scaling("cpu", motif_name)))
    return num_threads


def get_mem_mb_hmmer(motif_name):
    """Empirically, HMMER uses about 3g
    per CPU thread on average
    """

    num_threads = get_num_threads_hmmer(motif_name)
    mem_scale_factor = hmmer_scaling("mem", motif_name)

    mem_per_thread = 3072 * mem_scale_factor
    total_mem = int(mem_per_thread * num_threads)
    return total_mem


def hmmer_scaling(resource, motif_name):

    assert resource in ["cpu", "mem", "time"]
    try:
        scaling_factor = int(HMMER_MOTIF_SEARCH[motif_name][f"scale_{resource}"])
    except KeyError:
        scaling_factor = 1
    return scaling_factor


def hmmer_threshold_value(threshold, motif_name):

    # defaults as per HMMER cli interface
    _DEFAULT_HMMER_EVALUE = "10"
    _DEFAULT_HMMER_SCORE = 0

    assert threshold in ["evalue_t", "evalue", "score", "score_t"]

    t_value = None
    if threshold in ["evalue", "evalue_t"]:
        # NB: only the E-value threshold is used in the HMMER call
        t_value = HMMER_MOTIF_SEARCH[motif_name].get("evalue_t", _DEFAULT_HMMER_EVALUE)

    if threshold in ["score", "score_t"]:
        t_value = HMMER_MOTIF_SEARCH[motif_name].get("score_t", _DEFAULT_HMMER_EVALUE)

    assert t_value is not None

    return t_value
