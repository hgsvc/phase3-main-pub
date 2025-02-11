
rule asm_misjoins_paftools:
    """
    NB: rule could be augmented by adding a centromere
    annotation - of course, this does not necessarily
    exist in all case and has to be optional
    "paftools.js misjoin" parameter "-c FILE   BED for centromeres"
    """
    input:
        paf = rules.minimap_assembly_to_reference_align_paf.output.paf,
        cen_bed = lambda wildcards: load_reference_centromere_annotation(wildcards.refgenome)
    output:
        txt = DIR_PROC.joinpath(
            "60-flagging", "misjoins",
            "{sample}.{asm_unit}.{refgenome}.paftools-misjoins.txt"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    shell:
        "paftools.js misjoin -c {input.cen_bed} -e {input.paf} > {output.txt}"


rule run_all_paftools_misjoin:
    input:
        txt = expand(
            rules.asm_misjoins_paftools.output.txt,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN,
            refgenome=COMPLETE_REF_GENOME
        )
