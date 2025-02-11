
rule split_verkko_posthoc_phased_fasta:
    """
    This is a temporary rule to deal with
    the issue that Verkko does not split
    the main fasta into hap1/hap2 and
    unassigned if Verkko is restarted
    with the --paths argument.
    This should be fixed in some future
    release of Verkko.

    See gh#170
    """
    input:
        wait_file = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.wait"),
        path_check = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.paths.ok")
    output:
        check_file = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.ok")
    wildcard_constraints:
        phasing_state = "ps-sseq"
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt
    run:
        import io
        import os
        import pathlib as pl

        def is_empty(file_path):
            return os.stat(file_path).st_size == 0

        def dump_fasta_partition(part_name, part_seq_buffer, verkko_wd):
            part_file = verkko_wd.joinpath(f"assembly.{part_name}.fasta")
            assert part_file.is_file()
            assert is_empty(part_file)
            with open(part_file, "w") as fasta:
                _ = fasta.write(part_seq_buffer.getvalue())
            return

        verkko_run_wd = pl.Path(input.wait_file).with_suffix(".wd")
        assert verkko_run_wd.is_dir()

        main_assembly_file = verkko_run_wd.joinpath("assembly.fasta")
        assert main_assembly_file.is_file()
        assert not is_empty(main_assembly_file)

        expected_partitions = ["haplotype1", "haplotype2", "unassigned"]
        processed_partitions = set()
        active_partition = None
        partition_buffer = io.StringIO()
        num_contigs = 0
        with open(main_assembly_file, "r") as fasta:
            for line in fasta:
                if line.startswith(">"):
                    contig_partition = line.strip()[1:].split("-")[0]
                    num_contigs += 1
                    assert contig_partition in expected_partitions
                    if processed_partitions and contig_partition not in processed_partitions:
                        dump_fasta_partition(active_partition, partition_buffer, verkko_run_wd)
                        processed_partitions.add(active_partition)
                        partition_buffer = io.StringIO()
                        active_partition = contig_partition
                    active_partition = contig_partition
                    processed_partitions.add(active_partition)
                    partition_buffer.write(line)
                else:
                    partition_buffer.write(line)

        dump_fasta_partition(active_partition, partition_buffer, verkko_run_wd)
        processed_partitions.add(active_partition)
        assert sorted(processed_partitions) == expected_partitions, (
            f"Partitions: processed {processed_partitions} - expected {expected_partitions}"
        )

        with open(output.check_file, "w") as dump:
            _ = dump.write(f"Processed contigs: {num_contigs}\n")
    # END OF RUN BLOCK


localrules: collect_verkko_output_files
rule collect_verkko_output_files:
    input:
        check_file = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.ok")
    output:
        file_collection = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.output.json")
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script = find_script("collect_verkko_output"),
        verkko_wd = lambda wildcards, input: pathlib.Path(input.check_file).with_suffix(".wd"),
        base_dir = WORKDIR
    shell:
        "{params.script} --verkko-wd {params.verkko_wd} --base-dir {params.base_dir} "
            "--delete-logs --sample {wildcards.sample} --output {output.file_collection}"


################################
### 2023-08-03
# Verkko assembly output consists
# of duplicated sequences at least
# between disconnected and rDNA and
# disconnected and EBV. Unclear if
# discarding that would result in
# information loss, but for the
# time being, it is assumed that
# this is not the case. New rules
# added to deduplicate the output
# at least in the non-primary
# output files ("scraps").
###
#################################

rule compute_verkko_assembled_sequence_id:
    """
    The output computed here will eventually
    be discarded, and is only used to identify
    identical sequences among all Verkko output
    files (by sequence hash).
    """
    input:
        file_collection = DIR_PROC.joinpath("10-assemble/verkko/{sample}.{phasing_state}.output.json")
    output:
        stats = DIR_PROC.joinpath(
            "30-postprocess", "deduplicate", "seq_ids",
            "{sample}.{phasing_state}.verkko-asm-{asmtype}.statistics.tsv.gz"
        ),
    wildcard_constraints:
        asmtype = "(hap1|hap2|unassigned|disconnected|rdna|ebv|mito|wg)"
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        time_hrs=lambda wildcards, attempt: attempt,
    params:
        script=find_script("seqstats"),
        assembly=lambda wildcards, input: get_verkko_output(
            input.file_collection, f"{wildcards.asmtype}_fasta"
        ),
    shell:
        "{params.script} --cores {threads} "
        "--no-homopolymer-runs "
        "--no-short-tandem-repeats "
        "--no-sequence-composition "
        "--output-statistics {output.stats} "
        "--input-files {params.assembly}"


rule prepare_verkko_dup_seq_report:
    input:
        seq_stats = lambda wildcards: expand(
            rules.compute_verkko_assembled_sequence_id.output.stats,
            asmtype=get_verkko_asm_units(wildcards.phasing_state),
            allow_missing=True
        )
    output:
        report = DIR_RES.joinpath(
            "reports", "seq_dedup", "{sample}.{phasing_state}.verkko-dup-seq.tsv.gz"
        ),
        summary = DIR_RES.joinpath(
            "reports", "seq_dedup", "{sample}.{phasing_state}.verkko-dup-seq.summary.tsv"
        ),
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script = find_script("asm_dupseq_report")
    shell:
        "{params.script} --input {input.seq_stats} --report {output.report} "
            "--summary {output.summary}"


rule filter_verkko_dup_sequences:
    input:
        report = rules.prepare_verkko_dup_seq_report.output.report,
        file_collection = rules.collect_verkko_output_files.output.file_collection
    output:
        asm_unit = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.asm-{asm_unit}.fasta.gz"
        ),
        fai = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.asm-{asm_unit}.fasta.gz.fai"
        ),
        gzi = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.asm-{asm_unit}.fasta.gz.gzi"
        ),
    log:
        DIR_LOG.joinpath(
            "30-postprocess", "verkko",
            "{sample}.{phasing_state}.asm-{asm_unit}.dedup.log"
        ),
    wildcard_constraints:
        asm_unit = "(hap1|hap2|unassigned|disconnected|rdna|ebv|mito|wg)"
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script = find_script("filter_verkko_dupseq"),
        fasta = lambda wildcards, input: get_verkko_output(
            input.file_collection, f"{wildcards.asm_unit}_fasta"
        ),
        acc_res=lambda wildcards, output: register_result(output)
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "rm -f {output.asm_unit}.EMPTY ; "
        "if [ -s {params.fasta} ] ; then "
        "{{ "
        "{params.script} --input {params.fasta} --report {input.report} --verbose 2> {log}"
            " | "
        "bgzip > {output.asm_unit}"
            " && "
        "samtools faidx {output.asm_unit} ; "
        "}} else {{ "
        "touch {output.asm_unit} && touch {output.fai} "
        "&& touch {output.gzi} && touch {output.asm_unit}.EMPTY ; "
        "}} fi ;"


localrules: copy_verkko_exemplar_sequences
rule copy_verkko_exemplar_sequences:
    input:
        file_collection = rules.collect_verkko_output_files.output.file_collection
    output:
        ex_seq = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.exemplar-{asm_unit}.fasta.gz"
        ),
        fai = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.exemplar-{asm_unit}.fasta.gz.fai"
        ),
        gzi = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.exemplar-{asm_unit}.fasta.gz.gzi"
        ),
    wildcard_constraints:
        asm_unit = "(mito|ebv|rdna)"
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        fasta = lambda wildcards, input: get_verkko_output(
            input.file_collection, f"{wildcards.asm_unit}_repr"
        ),
        acc_res=lambda wildcards, output: register_result(output)
    shell:
        "rm -f {output.ex_seq}.EMPTY ; "
        "if [ -s {params.fasta} ] ; then "
        "{{ "
        "cat {params.fasta} | bgzip > {output.ex_seq}"
            " && "
        "samtools faidx {output.ex_seq} ; "
        "}} else {{ "
        "touch {output.ex_seq} && touch {output.fai} "
        "&& touch {output.gzi} touch {output.ex_seq}.EMPTY ; "
        "}} fi;"


### The following post-processing rules are an explicit code duplication
# because of the varying output that Verkko produces depending on
# the phasing mode. Maybe there will be more consistency in future versions
# of verkko

rule finalize_verkko_unphased_samples:
    """
    Rule for UNPHASED / ps-none samples
    """
    input:
        file_collection = rules.collect_verkko_output_files.output.file_collection
    output:
        hifi_cov = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "aux",
            "{sample}.{phasing_state}.node-hifi-cov.tsv"
        ),
        ont_cov = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "aux",
            "{sample}.{phasing_state}.node-ont-cov.tsv"
        ),
        gfa = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "graph",
            "{sample}.{phasing_state}.homopolymer-compressed-graph.gfa.gz"
        ),
        gfa_noseq = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "graph",
            "{sample}.{phasing_state}.homopolymer-compressed-graph.noseq.gfa"
        ),
        layout = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "graph",
            "{sample}.{phasing_state}.homopolymer-compressed-graph.layout.gz"
        )
    wildcard_constraints:
        phasing_state = "ps-none"
    params:
        acc_res=lambda wildcards, output: register_result(output),
        hifi_cov=lambda wildcards, input: get_verkko_output(input.file_collection, "hifi_cov"),
        ont_cov=lambda wildcards, input: get_verkko_output(input.file_collection, "ont_cov"),
        gfa=lambda wildcards, input: get_verkko_output(input.file_collection, "wg_gfa_hpc"),
        gfa_noseq=lambda wildcards, input: get_verkko_output(input.file_collection, "wg_gfa_noseq"),
        layout=lambda wildcards, input: get_verkko_output(input.file_collection, "wg_layout"),
    shell:
        "cp {params.hifi_cov} {output.hifi_cov} && "
        "cp {params.ont_cov} {output.ont_cov} && "
        "cat {params.gfa} | gzip > {output.gfa} && "
        "cp {params.gfa_noseq} {output.gfa_noseq} && "
        "cat {params.layout} | gzip > {output.layout}"


rule finalize_verkko_sseq_samples:
    """
    Rule for STRANDSEQ / ps-sseq samples
    """
    input:
        file_collection = rules.collect_verkko_output_files.output.file_collection,
        paths = lambda wildcards: MAP_SAMPLE_TO_INPUT_FILES[wildcards.sample]["phasing_paths"],
    output:
        layout = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "graph",
            "{sample}.{phasing_state}.homopolymer-compressed-graph.layout.gz"
        ),
        paths = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "aux",
            "{sample}.{phasing_state}.graphasing-paths.tsv"
        ),
    wildcard_constraints:
        phasing_state = "ps-sseq"
    params:
        acc_res=lambda wildcards, output: register_result(output),
        layout=lambda wildcards, input: get_verkko_output(input.file_collection, "wg_layout"),
    shell:
        "cat {params.layout} | gzip > {output.layout} && "
        "cp {input.paths} {output.paths}"


rule finalize_verkko_trio_samples:
    """
    Rule for TRIO / ps-trio samples
    """
    input:
        file_collection = rules.collect_verkko_output_files.output.file_collection
    output:
        hifi_cov = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "aux",
            "{sample}.{phasing_state}.node-hifi-cov.tsv"
        ),
        ont_cov = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "aux",
            "{sample}.{phasing_state}.node-ont-cov.tsv"
        ),
        gfa = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "graph",
            "{sample}.{phasing_state}.homopolymer-compressed-graph.gfa.gz"
        ),
        gfa_noseq = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "graph",
            "{sample}.{phasing_state}.homopolymer-compressed-graph.noseq.gfa"
        ),
        layout = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "graph",
            "{sample}.{phasing_state}.homopolymer-compressed-graph.layout.gz"
        ),
        coloring = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "aux",
            "{sample}.{phasing_state}.node-coloring.tsv"
        ),
        paths = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "aux",
            "{sample}.{phasing_state}.rukki-phasing-paths.tsv"
        ),
    wildcard_constraints:
        phasing_state = "ps-trio"
    params:
        acc_res=lambda wildcards, output: register_result(output),
        hifi_cov=lambda wildcards, input: get_verkko_output(input.file_collection, "hifi_cov"),
        ont_cov=lambda wildcards, input: get_verkko_output(input.file_collection, "ont_cov"),
        gfa=lambda wildcards, input: get_verkko_output(input.file_collection, "wg_gfa_hpc"),
        gfa_noseq=lambda wildcards, input: get_verkko_output(input.file_collection, "wg_gfa_noseq"),
        layout=lambda wildcards, input: get_verkko_output(input.file_collection, "wg_layout"),
        coloring=lambda wildcards, input: get_verkko_output(input.file_collection, "node_coloring"),
        paths=lambda wildcards, input: get_verkko_output(input.file_collection, "rukki_paths"),
    shell:
        "cp {params.hifi_cov} {output.hifi_cov} && "
        "cp {params.ont_cov} {output.ont_cov} && "
        "cat {params.gfa} | gzip > {output.gfa} && "
        "cp {params.gfa_noseq} {output.gfa_noseq} && "
        "cat {params.layout} | gzip > {output.layout} && "
        "cp {params.coloring} {output.coloring} && "
        "cp {params.paths} {output.paths}"


rule finalize_verkko_hic_samples:
    """
    Rule for HiC / ps-hic samples
    """
    input:
        file_collection = rules.collect_verkko_output_files.output.file_collection
    output:
        gfa_noseq = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "graph",
            "{sample}.{phasing_state}.homopolymer-compressed-graph.noseq.gfa"
        ),
        layout = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "graph",
            "{sample}.{phasing_state}.homopolymer-compressed-graph.layout.gz"
        ),
        paths = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "aux",
            "{sample}.{phasing_state}.hic-phasing-paths.tsv"
        ),
    wildcard_constraints:
        phasing_state = "ps-hic"
    params:
        acc_res=lambda wildcards, output: register_result(output),
        gfa_noseq=lambda wildcards, input: get_verkko_output(input.file_collection, "wg_gfa_noseq"),
        layout=lambda wildcards, input: get_verkko_output(input.file_collection, "wg_layout"),
        paths=lambda wildcards, input: get_verkko_output(input.file_collection, "rukki_paths"),
    shell:
        "cp {params.gfa_noseq} {output.gfa_noseq} && "
        "cat {params.layout} | gzip > {output.layout} && "
        "cp {params.paths} {output.paths}"


localrules: save_verkko_run_config
rule save_verkko_run_config:
    input:
        file_collection = rules.collect_verkko_output_files.output.file_collection
    output:
        run_yaml = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "aux",
            "{sample}.{phasing_state}.verkko-run-config.yml"
        )
    params:
        source_file=lambda wildcards, input: pathlib.Path(
            str(input.file_collection).replace("output.json", "wd")
        ).joinpath("verkko.yml"),
        acc_res=lambda wildcards, output: register_result(output)
    shell:
        "cp {params.source_file} {output.run_yaml}"


rule postprocess_verkko_unphased_samples:
    input:
        exemplars = expand(
            rules.copy_verkko_exemplar_sequences.output.ex_seq,
            phasing_state=["ps-none"],
            asm_unit=["rdna", "ebv", "mito"],
            sample=UNPHASED_SAMPLES
        ),
        asm_units = expand(
            rules.filter_verkko_dup_sequences.output.asm_unit,
            phasing_state=["ps-none"],
            asm_unit=get_verkko_asm_units("ps-none"),
            sample=UNPHASED_SAMPLES
        ),
        finalize = expand(
            rules.finalize_verkko_unphased_samples.output,
            sample=UNPHASED_SAMPLES,
            phasing_state=["ps-none"]
        ),
        runcfg = expand(
            rules.save_verkko_run_config.output.run_yaml,
            sample=UNPHASED_SAMPLES,
            phasing_state=["ps-none"]
        )


rule postprocess_verkko_sseq_samples:
    input:
        exemplars = expand(
            rules.copy_verkko_exemplar_sequences.output.ex_seq,
            phasing_state=["ps-sseq"],
            asm_unit=["rdna", "ebv", "mito"],
            sample=SSEQ_SAMPLES
        ),
        asm_units = expand(
            rules.filter_verkko_dup_sequences.output.asm_unit,
            phasing_state=["ps-sseq"],
            asm_unit=get_verkko_asm_units("ps-sseq"),
            sample=SSEQ_SAMPLES
        ),
        finalize = expand(
            rules.finalize_verkko_sseq_samples.output,
            sample=SSEQ_SAMPLES,
            phasing_state=["ps-sseq"]
        ),
        runcfg = expand(
            rules.save_verkko_run_config.output.run_yaml,
            sample=SSEQ_SAMPLES,
            phasing_state=["ps-sseq"]
        )


rule postprocess_verkko_trio_samples:
    input:
        exemplars = expand(
            rules.copy_verkko_exemplar_sequences.output.ex_seq,
            phasing_state=["ps-trio"],
            asm_unit=["rdna", "ebv", "mito"],
            sample=TRIO_SAMPLES
        ),
        asm_units = expand(
            rules.filter_verkko_dup_sequences.output.asm_unit,
            phasing_state=["ps-trio"],
            asm_unit=get_verkko_asm_units("ps-trio"),
            sample=TRIO_SAMPLES
        ),
        finalize = expand(
            rules.finalize_verkko_trio_samples.output,
            sample=TRIO_SAMPLES,
            phasing_state=["ps-trio"]
        ),
        runcfg = expand(
            rules.save_verkko_run_config.output.run_yaml,
            sample=TRIO_SAMPLES,
            phasing_state=["ps-trio"]
        )


rule postprocess_verkko_hic_samples:
    input:
        exemplars = expand(
            rules.copy_verkko_exemplar_sequences.output.ex_seq,
            phasing_state=["ps-hic"],
            asm_unit=["rdna", "ebv", "mito"],
            sample=HIC_SAMPLES
        ),
        asm_units = expand(
            rules.filter_verkko_dup_sequences.output.asm_unit,
            phasing_state=["ps-hic"],
            asm_unit=get_verkko_asm_units("ps-hic"),
            sample=HIC_SAMPLES
        ),
        finalize = expand(
            rules.finalize_verkko_hic_samples.output,
            sample=HIC_SAMPLES,
            phasing_state=["ps-hic"]
        ),
        runcfg = expand(
            rules.save_verkko_run_config.output.run_yaml,
            sample=HIC_SAMPLES,
            phasing_state=["ps-hic"]
        )
