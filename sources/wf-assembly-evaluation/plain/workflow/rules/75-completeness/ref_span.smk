
rule compute_approx_reference_span:
    input:
        aln_tsv = rules.normalize_mashmap_assembly_to_reference_align_paf.output.tsv,
        telo = rules.seqtk_annotate_telomere_motifs.output.table
    output:
        table = DIR_RES.joinpath(
            "coverage", "ref_span", "{refgenome}",
            "{sample}.{asm_unit}.{refgenome}.ref-span-telo.tsv"
        )
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        time_hrs=lambda wildcards, attempt: attempt
    params:
        script=find_script("ref_span.py")
    shell:
        "{params.script} --approx-align {input.aln_tsv} "
        "--telomere-bed {input.telo} --ref-span {output.table}"


rule run_all_compute_approx_ref_span:
    """TODO
    need an assembler-specific set of wildcards
    to the larger sequences / respective FASTAs
    """
    input:
        tables = expand(
            rules.compute_approx_reference_span.output.table,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN,
            refgenome=WILDCARDS_REF_GENOMES
        )
