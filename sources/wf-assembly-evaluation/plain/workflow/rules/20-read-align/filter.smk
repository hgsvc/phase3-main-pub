
rule filter_read_alignments_to_subset:
    """ This is designed to keep only
    the primary alignments (for tools
    such as NucFreq) or the supplementary
    alignments for downstream mosdepth
    windowed coverage computation.
    """
    input:
        bam = rules.align_reads_to_complete_assembly.output.mapped_bam,
        bai = rules.align_reads_to_complete_assembly.output.mapped_bai,
        clean_regions = get_clean_assembly_regions
        #clean_regions = rules.define_clean_assembly_regions.output.tag_tig
        # 2024-03-18: changed to realize skipping over contamination scan
    output:
        bam = DIR_PROC.joinpath(
            "20-read-align", "10_filter_per_file", "{sample}.{read_type}",
            "{sample}.{read_type}.{path_id}.{aln_subset}.sort.bam"
        ),
        bai = DIR_PROC.joinpath(
            "20-read-align", "10_filter_per_file", "{sample}.{read_type}",
            "{sample}.{read_type}.{path_id}.{aln_subset}.sort.bam.bai"
        ),
    benchmark:
        DIR_RSRC.joinpath(
            "20-read-align", "10_filter_per_file", "{sample}.{read_type}",
            "{sample}.{read_type}.{path_id}.{aln_subset}.samtools.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("aligner", "minimap.yaml")
    threads: CPU_LOW
    resources:
        mem_mb = lambda wildcards, attempt: {
            "onlySPL": 12288 + 12288 * attempt,
            "onlySEC": 12288 + 12288 * attempt,
            "onlyPRI": 24576 + 24576 * attempt
        }[wildcards.aln_subset],
        time_hrs = lambda wildcards, attempt: attempt * attempt,
        sort_mem_mb = lambda wildcards, attempt: 2048 + 2048 * attempt
    params:
        select_flag = lambda wildcards: {
            "onlyPRI": "--exclude-flags 3844",
            "onlySPL": "--require-flags 2048",
            "onlySEC": "--require-flags 256",
        }[wildcards.aln_subset],
        temp_prefix = lambda wildcards: f"temp/20-read-align/{wildcards.path_id}"
        # onlyPRI: unmapped; not primary; QC fail; duplicate; supplementary
        # onlySPL: supplementary
        # onlySEC: not primary
    shell:
        "mkdir -p {params.temp_prefix}"
            " && "
        "samtools view --threads {threads} -u {params.select_flag} "
        "--region-file {input.clean_regions} {input.bam} "
            " | "
        "samtools sort -l 9 -m {resources.sort_mem_mb}M --threads {threads} "
        "-T {params.temp_prefix}/{wildcards.sample}_{wildcards.read_type}_{wildcards.aln_subset}_mm2 "
        "-o {output.bam}"
            " && "
        "samtools index -@ {threads} {output.bam}"
            " ; "
        "rm -rfd {params.temp_prefix}"

