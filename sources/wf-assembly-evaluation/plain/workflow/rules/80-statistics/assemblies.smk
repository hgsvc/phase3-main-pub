

rule compute_assembly_sequence_statistics:
    input:
        fagz = get_asm_unit
    output:
        stats = DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.{asm_unit}.statistics.tsv.gz"
        ),
        summary = DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.{asm_unit}.summary.tsv"
        )
    benchmark:
        DIR_RSRC.joinpath("statistics", "assemblies", "{sample}.{asm_unit}.seqstats.rsrc")
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_HIGH
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt,
        time_hrs=lambda wildcards, attempt: attempt * attempt,
    params:
        script=find_script("seqstats"),
        report_seq_lens=lambda wildcards: SEQUENCE_LENGTH_THRESHOLDS_ASSEMBLY.get(wildcards.asm_unit, "default"),
        acc_res=lambda wildcards, output: register_result(output)
    shell:
        "{params.script} --cores {threads} "
        "--summary-length-thresholds {params.report_seq_lens} "
        "--temp-records 100 "
        "--str-motif-lengths 2 3 "
        "--output-statistics {output.stats} "
        "--output-summary {output.summary} "
        "--input-files {input.fagz}"


# TODO need a general solution for the assembly units
# that were specified in the sample sheet
rule run_all_assembly_statistics:
    input:
        stats = expand(
            rules.compute_assembly_sequence_statistics.output.summary,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_PLUS_CONTAM
        )
