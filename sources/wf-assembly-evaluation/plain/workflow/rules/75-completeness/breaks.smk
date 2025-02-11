
rule find_contig_alignment_breaks:
    """The output in query coordinates can be
    used as is in terms of an annotation file
    (dump to results folder).
    The output in reference coordinates must (should)
    be postprocessed to merge the results for all
    samples to then create a global annotation for
    remaining gaps in the assemblies.
    """
    input:
        norm_paf = get_contig_to_reference_norm_paf
    output:
        trg_all = DIR_PROC.joinpath(
            "75-completeness", "breaks", "{refgenome}",
            "{sample}.{asm_unit}.{refgenome}.{aln_type}.ctgaln-label.trg-all.bed.gz"
        ),
        trg_cmp = DIR_PROC.joinpath(
            "75-completeness", "breaks", "{refgenome}",
            "{sample}.{asm_unit}.{refgenome}.{aln_type}.ctgaln-label.trg-breaks.bed.gz"
        ),
        qry_all = DIR_PROC.joinpath(
            "75-completeness", "breaks", "{refgenome}",
            "{sample}.{asm_unit}.{refgenome}.{aln_type}.ctgaln-label.qry-all.bed.gz"
        ),
        qry_cmp = DIR_PROC.joinpath(
            "75-completeness", "breaks", "{refgenome}",
            "{sample}.{asm_unit}.{refgenome}.{aln_type}.ctgaln-label.qry-breaks.bed.gz"
        ),
    wildcard_constraints:
        aln_type = "(blevel|coarse)"
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("find_discontigs"),
        query_label=lambda wildcards: f"{wildcards.sample}.{wildcards.asm_unit}"
    shell:
        "{params.script} --alignments {input.norm_paf} --dump-bed-like "
        "--target-label {wildcards.refgenome} --query-label {params.query_label} "
        "--out-target-all {output.trg_all} --out-target-complement {output.trg_cmp} "
            " && "
        "{params.script} --alignments {input.norm_paf} --dump-bed-like "
        "--target-label {wildcards.refgenome} --query-label {params.query_label} "
        "--add-label-description "
        "--out-query-all {output.qry_all} --out-query-complement {output.qry_cmp} "


rule intersect_blevel_coarse_alignment_blocks:
    """Since all alignment blocks in the PAF files
    are identifiable via the MD5 hash, the orientation
    (target vs query) only marginally matters (most likely
    for the sequence/chromosome ends). We compute open/closed
    regions nevertheless for both views to increase confidence
    in final labeling.
    """
    input:
        blevel = DIR_PROC.joinpath(
            "75-completeness", "breaks", "{refgenome}",
            "{sample}.{asm_unit}.{refgenome}.blevel.ctgaln-label.{view}-all.bed.gz"
        ),
        coarse = DIR_PROC.joinpath(
            "75-completeness", "breaks", "{refgenome}",
            "{sample}.{asm_unit}.{refgenome}.coarse.ctgaln-label.{view}-all.bed.gz"
        ),
    output:
        isect = DIR_PROC.joinpath(
            "75-completeness", "breaks", "{refgenome}", "intersects",
            "{sample}.{asm_unit}.{refgenome}.ctgaln-label.{view}-all.isect.tsv.gz"
        ),
    wildcard_constraints:
        view="(trg|qry)"
    conda:
        DIR_ENVS.joinpath("biotools", "bedtools.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    shell:
        "bedtools intersect -wo -a {input.blevel} -b {input.coarse} | gzip > {output.isect}"


rule merge_gaps_alignment_blocks:
    input:
        qry_view = expand(
            rules.intersect_blevel_coarse_alignment_blocks.output.isect,
            view="qry",
            allow_missing=True
        ),
        qry_ngap = rules.merge_ngaps_annotations.output.bed,
        trg_view = expand(
            rules.intersect_blevel_coarse_alignment_blocks.output.isect,
            view="trg",
            allow_missing=True
        ),
        trg_ngap = rules.annotate_ngaps_in_reference.output.bed
    output:
        qry_view = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.{asm_unit}.{refgenome}.ctg-aln-gap.asm-coord.bed"
        ),
        trg_view = DIR_RES.joinpath(
            "regions", "{sample}",
            "{sample}.{asm_unit}.{refgenome}.ctg-aln-gap.ref-coord.bed"
        ),
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    params:
        script=find_script("merge_alignments")
    shell:
        "{params.script} --alignment-intersection {input.qry_view} --known-ngaps {input.qry_ngap} "
        "--ngap-label NGAPASM --add-label-description --output {output.qry_view}"
            " && "
        "{params.script} --alignment-intersection {input.trg_view} --known-ngaps {input.trg_ngap} "
        "--ngap-label NGAPREF --add-label-description --output {output.trg_view}"


rule run_all_label_contig_alignments:
    input:
        ctgaln_labeled = expand(
            rules.find_contig_alignment_breaks.output,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN,
            refgenome=WILDCARDS_REF_GENOMES,
            aln_type=["blevel", "coarse"]
        ),
        isect = expand(
            rules.intersect_blevel_coarse_alignment_blocks.output.isect,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN,
            refgenome=WILDCARDS_REF_GENOMES,
            view=["trg", "qry"]
        ),
        merged_gaps = expand(
            rules.merge_gaps_alignment_blocks.output,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN,
            refgenome=WILDCARDS_REF_GENOMES
        )
