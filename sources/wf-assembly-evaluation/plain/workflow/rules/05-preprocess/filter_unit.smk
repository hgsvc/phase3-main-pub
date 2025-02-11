
rule filter_asm_unit:
    """This rule was introduced to realize the scenario
    that the input assemblies are not checked for contamination
    but yet contain sequences that should be ignored during
    the evaluation.
    In hindsight, it would have been better to pre-filter
    all input assemblies for sequences to be hard-skipped
    for whatever reason
    TODO this is a candidate for code removal!
    NB: for samples going through this code path, 'skip_seqs'
    will be non-empty
    """
    input:
        fasta = lambda wildcards: SAMPLE_INFOS[wildcards.sample][
            ("asm", wildcards.asm_unit.split("-")[1], None)
        ],
        skip = lambda wildcards: SAMPLE_INFOS[wildcards.sample]["skip_seqs"]
    output:
        fasta = temp(
            DIR_PROC.joinpath(
                "05-preprocess", "filter_unit", "{sample}.{asm_unit}.filtered.fasta"
        ))
    log:
        DIR_LOG.joinpath("05-preprocess", "filter_unit", "{sample}.{asm_unit}.filtered.log")
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script = find_script("fasta_tag_merge"),
    shell:
        "{params.script} --input {input.fasta} --skip {input.skip} "
            "--report 2> {log} --output {output.fasta}"


rule compress_index_filtered_assembly_unit:
    input:
        fasta = rules.filter_asm_unit.output.fasta
    output:
        fagz = DIR_PROC.joinpath(
            "05-preprocess", "filter_unit", "{sample}.{asm_unit}.filtered.fa.gz"
        ),
        fai = DIR_PROC.joinpath(
            "05-preprocess", "filter_unit", "{sample}.{asm_unit}.filtered.fa.gz.fai"
        ),
        gzi = DIR_PROC.joinpath(
            "05-preprocess", "filter_unit", "{sample}.{asm_unit}.filtered.fa.gz.gzi"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "bgzip -c -l 6 -@ {threads} {input.fasta} > {output.fagz}"
            " && "
        "samtools faidx {output.fagz}"
