
rule homopolymer_compress_verkko_whole_genome:
    """NB: cannot use plain gzip/gunzip/zcat/pigz
    here because there is no built-in option to
    skip over/ignore empty files, which may happen
    for some files such as the mito sequences
    identified by Verkko (if executed w/ --screen).
    The usual "trick" of forcing a 0 return ('|| true')
    is dangerous because otherwise corrupted files
    are not detected here.
    TODO: seems like cleanest option is to run seq_hpc
    script once per input and to concatenate the coordinate
    maps afterwards. Needs a rewrite ...
    """
    input:
        fasta = lambda wildcards: expand(
            rules.filter_verkko_dup_sequences.output.asm_unit,
            asm_unit=get_verkko_asm_units(wildcards.phasing_state),
            allow_missing=True
        ),
        faidx = lambda wildcards: expand(
            rules.filter_verkko_dup_sequences.output.fai,
            asm_unit=get_verkko_asm_units(wildcards.phasing_state),
            allow_missing=True
        )
    output:
        fasta = DIR_PROC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.fastaseq.hpc.fasta.gz"
        ),
        faidx = DIR_PROC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.fastaseq.hpc.fasta.gz.fai"
        ),
        gzidx = DIR_PROC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.fastaseq.hpc.fasta.gz.gzi"
        ),
        cmap = DIR_PROC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.cmap.hpc.tsv.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.cmap.hpc.tsv.gz.tbi"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.fastaseq.hpc.rsrc"
        ),
    log:
        DIR_LOG.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.fastaseq.hpc.log"
        ),
    conda: DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt: 16384 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script = find_script("seq_hpc"),
        plain_tsv = lambda wildcards, output: pathlib.Path(output.cmap).with_suffix("")
    shell:
        "(pigz -c -d {input.fasta} || true) | {params.script} --cmap-table {params.plain_tsv} --report "
        "| bgzip -c > {output.fasta}"
            " && "
        "samtools faidx {output.fasta}"
            " && "
        "bgzip --threads {threads} {params.plain_tsv}"
            " && "
        "tabix --zero-based --sequence 1 --begin 2 --end 3 --comment \"#\" {output.cmap}"
            " ; "
        "rm -f {params.plain_tsv}"  # in case of rule failure


rule swap_hpc_to_plain_cmap:
    input:
        cmap = rules.homopolymer_compress_verkko_whole_genome.output.cmap
    output:
        cmap = DIR_PROC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.cmap.plain.tsv.gz"
        ),
        tbi = DIR_PROC.joinpath(
            "40-supplement", "verkko", "fasta_seq",
            "{sample}.{phasing_state}.cmap.plain.tsv.gz.tbi"
        )
    conda: DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt
    params:
        script = find_script("swap_cmap_orientation"),
    shell:
        "zcat {input.cmap} | {params.script} | bgzip -c --threads {threads} > {output.cmap}"
            " && "
        "tabix --zero-based --sequence 1 --begin 2 --end 3 --comment \"#\" {output.cmap}"


rule minimap_align_verkko_graphseq_to_fastaseq:
    input:
        fastaseq = rules.homopolymer_compress_verkko_whole_genome.output.fasta,
        gfaseq = get_verkko_gfaseq_hpc_fasta
    output:
        paf = DIR_PROC.joinpath(
            "40-supplement", "verkko", "gfa_to_fasta_align",
            "{sample}.{phasing_state}.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
    shell:
        "minimap2 -c -x asm5 -N 1 -p 0.95 --cs --eqx -t {threads} {input.fastaseq} {input.gfaseq}"
            " | "
        "pigz -p {threads} > {output.paf}"


rule normalize_minimap_gfa_to_fasta_align_paf:
    input:
        paf = rules.minimap_align_verkko_graphseq_to_fastaseq.output.paf
    output:
        tsv = DIR_PROC.joinpath(
            "40-supplement", "verkko", "gfa_to_fasta_align",
            "{sample}.{phasing_state}.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("normalize_paf")
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule expand_sequence_coordinates:
    input:
        paf = rules.normalize_minimap_gfa_to_fasta_align_paf.output.tsv,
        cmap = rules.homopolymer_compress_verkko_whole_genome.output.cmap,
        tbi = rules.homopolymer_compress_verkko_whole_genome.output.tbi
    output:
        tsv = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "aux",
            "{sample}.{phasing_state}.graph-linear-hpc-map.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt
    params:
        script=find_script("expand_coord"),
        acc_res=lambda wildcards, output: register_result(output)
    shell:
        "{params.script} --norm-paf {input.paf} --paf-coord-space hpc "
            "--coord-map {input.cmap} --paf-coord-expand target "
            "--out-table {output.tsv}"


rule run_verkko_unphased_supplement_cmap:
    input:
        tsv = expand(
            rules.expand_sequence_coordinates.output.tsv,
            sample=UNPHASED_SAMPLES,
            phasing_state=["ps-none"]
        )


rule run_verkko_sseq_supplement_cmap:
    """NB: this can only work if the unphased
    assembly run folder is still present because
    the above rule
    minimap_align_verkko_graphseq_to_fastaseq
    is desigend (see pyutils module) to select
    the assembly graph of the unphased run in
    that case.
    """
    input:
        tsv = expand(
            rules.expand_sequence_coordinates.output.tsv,
            sample=SSEQ_SAMPLES,
            phasing_state=["ps-sseq"]
        )


rule run_verkko_trio_supplement_cmap:
    input:
        tsv = expand(
            rules.expand_sequence_coordinates.output.tsv,
            sample=TRIO_SAMPLES,
            phasing_state=["ps-trio"]
        )


if False:
    # NB: Verkko does not output the complete
    # assembly graph including sequences in these
    # two cases, hence, no coordinate mapping can
    # be produced. Note that for Strand-seq, the
    # respective coordinate mapping is available
    # via the unphased assembly run.
    rule run_verkko_hic_supplement_cmap:
        input:
            tsv = expand(
                rules.expand_sequence_coordinates.output.tsv,
                sample=HIC_SAMPLES,
                phasing_state=["ps-hic"]
            )
