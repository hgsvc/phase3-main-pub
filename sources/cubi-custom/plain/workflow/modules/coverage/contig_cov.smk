
rule verkko_dump_assembly_coverage_issue_table:
    input:
        tsv = WORKDIR_EVAL.joinpath(
            "results/coverage/contig_ref/t2tv2",
            "{sample}.t2tv2.win-ctg-cov.tsv.gz"
        ),
        aln_hap1 = WORKDIR_EVAL.joinpath(
            "results/alignments/contig_to_ref/t2tv2/table",
            "{sample}.asm-hap1.t2tv2.norm-paf.tsv.gz"
        ),
        aln_hap2 = WORKDIR_EVAL.joinpath(
            "results/alignments/contig_to_ref/t2tv2/table",
            "{sample}.asm-hap2.t2tv2.norm-paf.tsv.gz"
        ),
        aln_unassigned = WORKDIR_EVAL.joinpath(
            "results/alignments/contig_to_ref/t2tv2/table",
            "{sample}.asm-unassigned.t2tv2.norm-paf.tsv.gz"
        ),
        aln_disconnected = WORKDIR_EVAL.joinpath(
            "results/alignments/contig_to_ref/t2tv2/table",
            "{sample}.asm-disconnected.t2tv2.norm-paf.tsv.gz"
        ),
    output:
        tview_bed = DIR_RES.joinpath(
            "issues", "contig_ref", "t2tv2",
            "{sample}.t2tv2.ctg-cov-issues.ref.bed.gz"
        ),
        qview_bed = DIR_RES.joinpath(
            "issues", "contig_ref", "t2tv2",
            "{sample}.t2tv2.ctg-cov-issues.assm.bed.gz"
        ),
        stats = DIR_RES.joinpath(
            "issues", "contig_ref", "t2tv2",
            "{sample}.t2tv2.ctg-cov-issues.stats.tsv"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script = DIR_SCRIPTS.joinpath("issue_tracks.py"),
        tmp_tview = lambda wildcards, output: pathlib.Path(output.tview_bed).with_suffix(".tmp.bed"),
        tmp_qview = lambda wildcards, output: pathlib.Path(output.qview_bed).with_suffix(".tmp.bed"),
        acros = ACROS
    shell:
        "{params.script} --contig-cov {input.tsv} --aln-unassign {input.aln_unassigned} "
        "--aln-hap1 {input.aln_hap1} --aln-hap2 {input.aln_hap2} --aln-disconnect {input.aln_disconnected} "
        "--acrocentrics {params.acros} --window-size 10000 --out-stats {output.stats} "
        "--out-target-view {params.tmp_tview} --out-query-view {params.tmp_qview}"
            " && "
        "bgzip -c {params.tmp_tview} > {output.tview_bed}"
            " && "
        "tabix -p bed {output.tview_bed} && rm {params.tmp_tview}"
            " && "
        "bgzip -c {params.tmp_qview} > {output.qview_bed}"
            " && "
        "tabix -p bed {output.qview_bed} && rm {params.tmp_qview}"


rule run_all_issue_tracks:
    input:
        stats = expand(
            rules.verkko_dump_assembly_coverage_issue_table.output.stats,
            sample=SAMPLES
        )
