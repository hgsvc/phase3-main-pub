
rule seqtk_annotate_telomere_motifs:
    input:
        fasta = get_asm_unit
    output:
        table = DIR_PROC.joinpath(
            "70-annotate", "telomeres", "seqtk",
            "{sample}.{asm_unit}.telo.tsv"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        time_hrs=lambda wildcards, attempt: attempt
    conda:
        DIR_ENVS.joinpath("biotools", "seqtk.yaml")
    shell:
        "seqtk telo {input.fasta} > {output.table}"


rule run_all_seqtk_telomere_motifs:
    """TODO introduce assembler-specific
    parameter for asm_unit wildcards
    """
    input:
        tables = expand(
            rules.seqtk_annotate_telomere_motifs.output.table,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN
        )

