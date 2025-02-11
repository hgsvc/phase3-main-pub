
rule dump_verkko_gfa_sequences:
    """
    """
    input:
        file_collection = rules.collect_verkko_output_files.output.file_collection
    output:
        fasta = DIR_PROC.joinpath(
            "40-supplement", "verkko", "graph_seq",
            "{sample}.{phasing_state}.gfaseq.hpc.fasta.gz"
        ),
        faidx = DIR_PROC.joinpath(
            "40-supplement", "verkko", "graph_seq",
            "{sample}.{phasing_state}.gfaseq.hpc.fasta.gz.fai"
        ),
    conda:
        DIR_ENVS.joinpath("graphtools.yaml")
    wildcard_constraints:
        phasing_state = "(ps-none|ps-trio)"
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt
    params:
        gfa = lambda wildcards, input: get_verkko_output(input.file_collection, "wg_gfa_hpc")
    shell:
        "gfatools gfa2fa -l 0 {params.gfa} | bgzip -c > {output.fasta}"
            " && "
        "samtools faidx {output.fasta}"
