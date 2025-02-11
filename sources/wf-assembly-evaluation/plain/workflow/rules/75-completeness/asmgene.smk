
rule ref_completeness_genemodel_any:
    """
    This rule and the two subsequent ones violate the principle
    of not duplicating code. Why? Creating sex-specific reference
    gene model alignments in a single rule would require intro-
    ducing a new wildcard (sex, e.g., any / male / female) but the
    karyotype assignment for phased assemblies happens as part of
    this workflow. Hence, the information which sex-specific gene
    model to use is only available as an output of this workflow.
    Feeding back that information as a new wildcard to trigger
    all runs in the rule
    > run_all_asm_completeness_asmgene
    is not explicitly possible. As a consequence, duplicating the
    code for the reference gene model alignment seems to be the
    cleaner solution, which also enables a straightforward selection
    of the right gene model via an input function in the rule
    > asm_completeness_asmgene_stats
    """
    input:
        fasta = lambda wildcards: get_reference_genome(wildcards.refgenome),
        cdna = lambda wildcards: get_gene_model(wildcards.genemodel, None, wildcards.refgenome, "any")
    output:
        paf = DIR_PROC.joinpath(
            "75-completeness", "asmgene",
            "{refgenome}.{genemodel}.any.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
    shell:
        "minimap2 -cxsplice:hq -t {threads} {input.fasta} {input.cdna}"
            " | "
        "pigz > {output.paf}"


rule ref_completeness_genemodel_female:
    input:
        fasta = lambda wildcards: get_reference_genome(wildcards.refgenome, "female-hap"),
        cdna = lambda wildcards: get_gene_model(wildcards.genemodel, None, wildcards.refgenome, "female")
    output:
        paf = DIR_PROC.joinpath(
            "75-completeness", "asmgene",
            "{refgenome}.{genemodel}.female.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
    shell:
        "minimap2 -cxsplice:hq -t {threads} {input.fasta} {input.cdna}"
            " | "
        "pigz > {output.paf}"


rule ref_completeness_genemodel_male:
    input:
        fasta = lambda wildcards: get_reference_genome(wildcards.refgenome, "male-hap"),
        cdna = lambda wildcards: get_gene_model(wildcards.genemodel, None, wildcards.refgenome, "male")
    output:
        paf = DIR_PROC.joinpath(
            "75-completeness", "asmgene",
            "{refgenome}.{genemodel}.male.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
    shell:
        "minimap2 -cxsplice:hq -t {threads} {input.fasta} {input.cdna}"
            " | "
        "pigz > {output.paf}"


rule asm_completeness_genemodel:
    """
    NB: the global variable COMPLETE_REF_GENOME is used
    here to avoid adding a new wildcard ("refgenome") that
    is not directly related to the gene model to assembly
    alignment.
    """
    input:
        asm_karyo = expand(
            rules.estimate_asm_unit_karyotype.output.karyo_est,
            refgenome=COMPLETE_REF_GENOME,
            allow_missing=True
        ),
        fasta = get_asm_unit,
        cdna = lambda wildcards: get_gene_model(wildcards.genemodel, wildcards.sample, COMPLETE_REF_GENOME, wildcards.asm_unit)
    output:
        paf = DIR_PROC.joinpath(
            "75-completeness", "asmgene",
            "assemblies", "{sample}",
            "{sample}.{asm_unit}.{genemodel}.paf.gz"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "75-completeness", "asmgene",
            "assemblies",
            "{sample}.{asm_unit}.{genemodel}.mm2.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
    shell:
        "minimap2 -cxsplice:hq -t {threads} {input.fasta} {input.cdna}"
            " | "
        "pigz > {output.paf}"


rule asm_completeness_asmgene_stats:
    """
    paftools.js asmgene -e == print fragmented/missing genes
    """
    input:
        asm_karyo = rules.estimate_asm_unit_karyotype.output.karyo_est,
        ref = lambda wildcards: get_reference_gene_model_alignment(
            wildcards.genemodel, wildcards.sample, wildcards.refgenome, wildcards.asm_unit
        ),
        assm = rules.asm_completeness_genemodel.output.paf
    output:
        stats = DIR_PROC.joinpath(
            "75-completeness", "asmgene",
            "assemblies", "{sample}",
            "{sample}.{asm_unit}.{refgenome}.{genemodel}.asmgene-stats.txt"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: 1
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt - 1,
        mem_mb = lambda wildcards, attempt: 2048 * attempt,
    shell:
        "paftools.js asmgene -e {input.ref} {input.assm} > {output.stats}"


rule asm_completeness_asmgene_process_output:
    input:
        txt = rules.asm_completeness_asmgene_stats.output.stats
    output:
        issues = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.{asm_unit}.{refgenome}.{genemodel}.asmgene-issues.bed"
        ),
        stats = DIR_RES.joinpath(
            "statistics", "completeness", "asmgene",
            "{sample}.{asm_unit}.{refgenome}.{genemodel}.asmgene-stats.tsv"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    params:
        script=find_script("split_asmgene_stats"),
        acc_res=lambda wildcards, output: register_result(output.issues, output.stats)
    shell:
        "{params.script} --input {input.txt} "
        "--bed-issues {output.issues} --stats-out {output.stats}"


rule run_all_asm_completeness_asmgene:
    input:
        asmgene_stats = expand(
            rules.asm_completeness_asmgene_stats.output.stats,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            refgenome=COMPLETE_REF_GENOME,
            genemodel=WILDCARDS_GENE_MODELS
        ),
        karyotypes = rules.build_assembly_karyotype_summary.output.tsv,
        agg_stats = expand(
            rules.asm_completeness_asmgene_process_output.output.stats,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN,
            refgenome=COMPLETE_REF_GENOME,
            genemodel=WILDCARDS_GENE_MODELS
        )
