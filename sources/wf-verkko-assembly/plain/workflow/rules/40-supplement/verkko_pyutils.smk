
def get_verkko_gfaseq_hpc_fasta(wildcards):
    """
    The Strand-seq phasing of Verkko graphs
    requires restarting Verkko, which implies
    that the corresponding GFA (assembly graph)
    is actually from the unphased assembly.
    This function exists to properly map between
    phasing state wildcard and the required
    graph sequences.
    """

    known_phasing_states = {
        "ps-sseq": DIR_PROC.joinpath(
            "40-supplement", "verkko", "graph_seq",
            "{sample}.ps-none.gfaseq.hpc.fasta.gz"
        ),
        "ps-trio": DIR_PROC.joinpath(
            "40-supplement", "verkko", "graph_seq",
            "{sample}.ps-trio.gfaseq.hpc.fasta.gz"
        ),
        "ps-none": DIR_PROC.joinpath(
            "40-supplement", "verkko", "graph_seq",
            "{sample}.ps-none.gfaseq.hpc.fasta.gz"
        ),
    }

    _this_func = "40-supplement::verkko_pyutils::get_verkko_gfaseq_hpc_fasta"

    if not hasattr(wildcards, "phasing_state"):
        logerr(f"{_this_func}: no wildcard 'phasing state' --- {wildcards}")
        raise ValueError(wildcards)
    if wildcards.phasing_state not in known_phasing_states:
        logerr(f"{_this_func}: unknown phasing state --- {wildcards.phasing_state}")
        raise ValueError(wildcards)

    return known_phasing_states[wildcards.phasing_state]
