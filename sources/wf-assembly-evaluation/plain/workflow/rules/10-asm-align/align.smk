
rule minimap_assembly_to_reference_align_paf:
    """NB: given the wildcard composition,
    the rules in this module are not supposed
    to be used to align sequence contaminants
    to the reference (which seems like an
    illogical thing to do anyway).
    """
    input:
        ref = lambda wildcards: get_reference_assembly(wildcards.sample, wildcards.refgenome),
        assm = get_asm_unit
    output:
        paf = DIR_RES.joinpath(
            "alignments", "contig_to_ref", "{refgenome}",
            "paf", "{sample}.{asm_unit}.{refgenome}.paf.gz"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
    shell:
        "minimap2 -c -x asm20 --cs --eqx -t {threads} {input.ref} {input.assm}"
            " | "
        "pigz -p {threads} > {output.paf}"


rule normalize_minimap_assembly_to_reference_align_paf:
    input:
        paf = rules.minimap_assembly_to_reference_align_paf.output.paf
    output:
        tsv = DIR_RES.joinpath(
            "alignments", "contig_to_ref", "{refgenome}",
            "table", "{sample}.{asm_unit}.{refgenome}.norm-paf.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("normalize_paf")
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


rule minimap_assembly_to_reference_align_bam:
    input:
        ref = lambda wildcards: get_reference_assembly(wildcards.sample, wildcards.refgenome),
        assm = get_asm_unit
    output:
        bam = DIR_RES.joinpath(
            "alignments", "contig_to_ref", "{refgenome}",
            "bam", "{sample}.{asm_unit}.{refgenome}.sort.bam"
        ),
        bai = DIR_RES.joinpath(
            "alignments", "contig_to_ref", "{refgenome}",
            "bam", "{sample}.{asm_unit}.{refgenome}.sort.bam.bai"
        ),
        unmapped = DIR_RES.joinpath(
            "alignments", "contig_to_ref", "{refgenome}",
            "bam", "{sample}.{asm_unit}.{refgenome}.unmapped.bam"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_MEDIUM
    resources:
        time_hrs = lambda wildcards, attempt: 1 * attempt,
        mem_mb = lambda wildcards, attempt: 32768 + 32768 * attempt,
        sort_mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        readgroup = lambda wildcards: (
            f'"@RG\\tID:{wildcards.sample}_{wildcards.asm_unit}'
            f'\\tSM:{wildcards.sample}"'
        ),
        sam_flag_out = 1540,  # unmap, PCR-dup, QC-fail --- keep 2nd align!
        sam_threads = CPU_LOW,
    shell:
        "minimap2 -a -x asm20 --cs --eqx -t {threads}"
        " -R {params.readgroup} {input.ref} {input.assm}"
            " | "
        " samtools view -u -h --output-unselected {output.unmapped} "
        " -F {params.sam_flag_out} --threads {params.sam_threads}"
            " | "
        " samtools sort -l 9 -m {resources.sort_mem_mb}M "
        " --threads {params.sam_threads} "
        " -T {wildcards.sample}_{wildcards.asm_unit}_{wildcards.refgenome}_mm2 -o {output.bam} "
            " && "
        "samtools index -@ {threads} {output.bam}"


rule mashmap_assembly_to_reference_align_paf:
    """The coarse-grained mashmap alignments are
    only used to provide a rough estimate of
    reference/chromosome spanning in particular
    for incomplete reference genomes (i.e. hg38)

    TODO - make cli parameters user-configurable
    --pi 99 / percent seq. identity
    --segLength 100000 / skip shorter sequences
    """
    input:
        ref = lambda wildcards: get_reference_assembly(wildcards.sample, wildcards.refgenome),
        assm = get_asm_unit
        #idx = rules.compress_clean_assembly_sequences.output.fai,
        # existence of index is by construction or is enforced
        # by assembly lookup function / 2024-03-18
    output:
        paf = DIR_PROC.joinpath(
            "10-asm-align", "contig_to_ref", "{refgenome}",
            "paf", "{sample}.{asm_unit}.{refgenome}.approx.paf"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "mashmap.yaml")
    threads: CPU_LOW
    resources:
        time_hrs = lambda wildcards, attempt: attempt,
        mem_mb = lambda wildcards, attempt: 32768 * attempt,
    shell:
        "mashmap -r {input.ref} -q {input.assm} "
        "-f one-to-one --pi 99 --segLength 100000 --dense "
        "--output {output.paf}"


rule normalize_mashmap_assembly_to_reference_align_paf:
    input:
        paf = rules.mashmap_assembly_to_reference_align_paf.output.paf
    output:
        tsv = DIR_PROC.joinpath(
            "10-asm-align", "contig_to_ref", "{refgenome}",
            "table", "{sample}.{asm_unit}.{refgenome}.norm-approx.tsv.gz"
        )
    conda:
        DIR_ENVS.joinpath("pyutils.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt
    params:
        script=find_script("normalize_paf")
    shell:
        "{params.script} --input {input.paf} --output {output.tsv}"


# TODO need general strategy to identify assembly units to use here
rule run_minimap_contig_to_ref_alignments:
    input:
        bams = expand(
                rules.minimap_assembly_to_reference_align_bam.output.bam,
                refgenome=WILDCARDS_REF_GENOMES,
                sample=SAMPLES,
                asm_unit=ASSEMBLY_UNITS_NO_CONTAM
        ),
        paf = expand(
                rules.normalize_minimap_assembly_to_reference_align_paf.output.tsv,
                refgenome=WILDCARDS_REF_GENOMES,
                sample=SAMPLES,
                asm_unit=ASSEMBLY_UNITS_NO_CONTAM
        )


rule run_mashmap_contig_to_ref_alignments:
    """TODO
    This needs assembler-specific wildcards;
    only larger contigs make sense for the approx
    alignment
    """
    input:
        tsv = expand(
            rules.normalize_mashmap_assembly_to_reference_align_paf.output.tsv,
            refgenome=WILDCARDS_REF_GENOMES,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_MAIN
        )
