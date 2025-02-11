
rule annotate_ngaps_in_reference:
    """This rule was introduced to ensure that also
    N gaps in the reference (hg38 etc.) are properly
    considered in the module
    75-completeness::breaks.smk
    """
    input:
        fasta = lambda wildcards: DIR_GLOBAL_REF.joinpath(config["refgenomes"][wildcards.refgenome]["any"])
    output:
        bed = DIR_LOCAL_REF.joinpath("{refgenome}.ngaps.bed")
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("localize_ngaps")
    shell:
        "{params.script} --fasta-input {input.fasta} "
        "--output {output.bed} --name {wildcards.refgenome}"


rule generate_ngaps_annotation:
    input:
        asm_unit = get_asm_unit
    output:
        bed = temp(
            DIR_PROC.joinpath(
                "05-preprocess", "ngaps", "{sample}",
                "{sample}.{asm_unit}.ngaps.bed"
            )
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("localize_ngaps")
    shell:
        "{params.script} --fasta-input {input.asm_unit} "
        "--output {output.bed} --name {wildcards.sample}"


localrules: merge_ngaps_annotations
rule merge_ngaps_annotations:
    input:
        tables = expand(
            rules.generate_ngaps_annotation.output.bed,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            allow_missing=True
        )
    output:
        bed = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.ngaps.bed"
        )
    run:
        import pandas as pd
        merged = []
        for table in input.tables:
            df = pd.read_csv(table, sep="\t", header=0)
            merged.append(df)
        merged = pd.concat(merged, axis=0, ignore_index=False)
        merged.sort_values(["#contig", "start"], inplace=True)
        merged.to_csv(output.bed, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule run_all_ngaps_annotation:
    input:
        beds = expand(
            rules.merge_ngaps_annotations.output.bed,
            sample=SAMPLES
        )
