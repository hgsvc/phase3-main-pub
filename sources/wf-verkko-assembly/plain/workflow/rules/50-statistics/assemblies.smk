
rule compute_verkko_assembly_stats:
    input:
        asm_unit = DIR_RES.joinpath(
            "assemblies", "verkko", "{sample}.{phasing_state}",
            "{sample}.{phasing_state}.asm-{asm_unit}.fasta.gz"
        ),
    output:
        stats = DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.{phasing_state}.verkko-asm-{asm_unit}.statistics.tsv.gz"
        ),
        summary = DIR_RES.joinpath(
            "statistics", "assemblies", "{sample}.{phasing_state}.verkko-asm-{asm_unit}.summary.tsv"
        )
    benchmark:
        DIR_RSRC.joinpath("statistics", "assemblies", "{sample}.{phasing_state}.verkko-asm-{asm_unit}.stats.rsrc")
    wildcard_constraints:
        asmtype = "(hap1|hap2|unassigned|rdna)"
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_HIGH
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt,
        time_hrs=lambda wildcards, attempt: attempt * attempt,
    params:
        script=find_script("seqstats"),
        report_seq_lens = " ".join(map(str, [int(1e5), int(5e5), int(1e6), int(1e7), int(5e7), int(1e8)])),
        acc_res=lambda wildcards, output: register_result(output.stats, output.summary)
    shell:
        "{params.script} --cores {threads} "
        "--summary-length-thresholds {params.report_seq_lens} "
        "--no-homopolymer-runs "
        "--no-canonical-sequence "
        "--temp-records 100 "
        "--str-motif-lengths 2 3 "
        "--output-statistics {output.stats} "
        "--output-summary {output.summary} "
        "--input-files {input.asm_unit}"


rule get_verkko_unphased_output_stats:
    input:
        summary = expand(
            DIR_RES.joinpath(
                "statistics", "assemblies",
                "{sample}.ps-none.verkko-{vrk_out}.summary.tsv"
            ),
            sample=UNPHASED_SAMPLES,
            vrk_out=["asm-wg", "asm-rdna"]
        )


rule get_verkko_sseq_phased_output_stats:
    input:
        summary = expand(
            DIR_RES.joinpath(
                "statistics", "assemblies",
                "{sample}.ps-sseq.verkko-{vrk_out}.summary.tsv"
            ),
            sample=SSEQ_SAMPLES,
            vrk_out=["asm-hap1", "asm-hap2", "asm-unassigned", "asm-rdna"]
        )


rule get_verkko_trio_phased_output_stats:
    input:
        summary = expand(
            DIR_RES.joinpath(
                "statistics", "assemblies",
                "{sample}.ps-trio.verkko-{vrk_out}.summary.tsv"
            ),
            sample=TRIO_SAMPLES,
            vrk_out=["asm-hap1", "asm-hap2", "asm-unassigned", "asm-rdna"]
        )


rule get_verkko_hic_phased_output_stats:
    input:
        summary = expand(
            DIR_RES.joinpath(
                "statistics", "assemblies",
                "{sample}.ps-hic.verkko-{vrk_out}.summary.tsv"
            ),
            sample=HIC_SAMPLES,
            vrk_out=["asm-hap1", "asm-hap2", "asm-unassigned", "asm-rdna"]
        )
