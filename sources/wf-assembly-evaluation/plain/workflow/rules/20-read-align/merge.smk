
rule merge_read_to_assembly_subset_alignments:
    input:
        bams = lambda wildcards: expand(
            rules.filter_read_alignments_to_subset.output.bam,
            path_id=SAMPLE_INFOS[wildcards.sample][("reads", wildcards.read_type, "path_id_all")],
            allow_missing=True
        ),
        bai = lambda wildcards: expand(
            rules.filter_read_alignments_to_subset.output.bai,
            path_id=SAMPLE_INFOS[wildcards.sample][("reads", wildcards.read_type, "path_id_all")],
            allow_missing=True
        )
    output:
        bam = DIR_PROC.joinpath(
            "20-read-align", "20_merge_subsets", "{sample}.{read_type}",
            "{sample}.{read_type}.{aln_subset}.sort.bam"
        ),
        bai = DIR_PROC.joinpath(
            "20-read-align", "20_merge_subsets", "{sample}.{read_type}",
            "{sample}.{read_type}.{aln_subset}.sort.bam.bai"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            "20-read-align", "20_merge_subsets", "{sample}.{read_type}",
            "{sample}.{read_type}.{aln_subset}.samtools.rsrc"
        )
    wildcard_constraints:
        aln_subset = "(" + "|".join(["onlyPRI", "onlySPL", "onlySEC"]) + ")"
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt: 1024 * attempt,
        time_hrs = lambda wildcards, attempt: attempt,
    shell:
        "samtools merge -r -f --threads {threads} -o {output.bam} {input.bams}"
            " && "
        "samtools index -@ {threads} {output.bam}"


rule extract_primary_alignment_read_lists:
    input:
        bam = rules.merge_read_to_assembly_subset_alignments.output.bam,
        bai = rules.merge_read_to_assembly_subset_alignments.output.bai
    output:
        read_list = DIR_RES.joinpath(
            "read_sets", "{sample}.{read_type}.{aln_subset}.reads.tsv.gz"
        ),
        cov_cache = DIR_PROC.joinpath(
            "20-read-align", "30_cache_cov", "{sample}.{read_type}.{aln_subset}.ctg-cov.h5"
        )
    conda:
        DIR_ENVS.joinpath("pyseq.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 4096 + 4096 * attempt,
        time_hrs=lambda wildcards, attempt: attempt
    params:
        script=find_script("get_pct_aligned"),
        min_length=lambda wildcards: {"hifi": 5000, "ont": 10000}[wildcards.read_type],
        min_aligned=lambda wildcards: {"hifi": 99, "ont": 95}[wildcards.read_type]
    shell:
        "{params.script} --input {input.bam} --io-threads {threads} "
            "--min-pct-aligned {params.min_aligned} --min-read-length {params.min_length} "
            "-tsv {output.read_list} -hdf {output.cov_cache}"


rule run_all_get_primary_alignment_read_lists:
    """
    TODO
    read types should be encoded as wildcard variable
    """
    input:
        tsv = expand(
            rules.extract_primary_alignment_read_lists.output.read_list,
            sample=SAMPLES,
            aln_subset=["onlyPRI"],
            read_type=["hifi", "ont"]
        )
