
rule prepare_region_cache:
    input:
        sizes = DIR_PROC.joinpath("asm_label_qc", "assembly_size", "{sample}.sizes.txt"),
        flagger = rules.normalize_flagger_annotation.output.bed_like,
        flagbin = rules.binarize_flagger_subset.output.bed,
        inspect_hifi = DIR_RES.joinpath("asm_label_qc", "norm_tables",
            "inspector", "{sample}.hifi.inspector-errors.tsv.gz"
        ),
        inspect_ont = DIR_RES.joinpath("asm_label_qc", "norm_tables",
            "inspector", "{sample}.ont.inspector-errors.tsv.gz"
        ),
        merqury = rules.normalize_merqury_kmer_tracks.output.bed_like,
        hifi_mq00 = DIR_RES.joinpath("asm_label_qc", "norm_tables", "read_depth",
            "{sample}.hifi.mq00.mosdepth-windowed.tsv.gz"
        ),
        hifi_mq60 = DIR_RES.joinpath("asm_label_qc", "norm_tables", "read_depth",
            "{sample}.hifi.mq60.mosdepth-windowed.tsv.gz"
        ),
        ont_mq00 = DIR_RES.joinpath("asm_label_qc", "norm_tables", "read_depth",
            "{sample}.ont.mq00.mosdepth-windowed.tsv.gz"
        ),
        ont_mq60 = DIR_RES.joinpath("asm_label_qc", "norm_tables", "read_depth",
            "{sample}.ont.mq60.mosdepth-windowed.tsv.gz"
        ),
        nucfreq = rules.bin_nucfreq_regions_by_coverage.output.bed,
        nfbin = rules.binarize_nucfreq_output.output.bed,
        sd95 = rules.split_segdup_annotation.output.sd95,
        sd98 = rules.split_segdup_annotation.output.sd98,
        cen = DIR_RES.joinpath(
            "asm_label_qc", "norm_tables", "centromeres",
            "{sample}.active_asat_HOR_arrays_v3.bed"
        ),
        sseqbrkp = rules.translate_hpc_breakpoint_coordinates.output.bed,
        busco = rules.normalize_busco_issue_annotation.output.bed
    output:
        hdf = DIR_RES.joinpath(
            "asm_label_qc", "region_cache",
            "{sample}.assessem-cache.h5"
        )
    log:
        DIR_LOG.joinpath(
            "asm_label_qc", "region_cache", "{sample}.assessem-cache.log"
        )
    conda:
        DIR_ENVS.joinpath("assessem.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt
    params:
        script=DIR_SCRIPTS.joinpath("asm_label_qc", "assessem-cache.py"),
        track_files = lambda wildcards, input: get_assessem_cli_parameters(input, "files"),
        track_labels = lambda wildcards, input: get_assessem_cli_parameters(input, "labels"),
        score_columns = lambda wildcards, input: get_assessem_cli_parameters(input, "columns"),
    shell:
        "{params.script} --verbose --data-cache {output.hdf} --genome-size {input.sizes} "
        " --track-files {params.track_files} --track-labels {params.track_labels} "
        " --score-columns {params.score_columns} > {log}"


rule compute_embedding:
    input:
        cache = rules.prepare_region_cache.output.hdf
    output:
        embed = DIR_RES.joinpath(
            "asm_label_qc", "embedding",
            "{sample}.{bin_size}.{feature_set}.embed-dataset.npy"
        ),
        dataset = DIR_RES.joinpath(
            "asm_label_qc", "embedding",
            "{sample}.{bin_size}.{feature_set}.binned-dataset.tsv.gz"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "asm_label_qc", "embedding", "{sample}.{bin_size}.{feature_set}.assessem-embed.rsrc"
        )
    log:
        DIR_LOG.joinpath(
            "asm_label_qc", "embedding", "{sample}.{bin_size}.{feature_set}.assessem-embed.log"
        )
    conda:
        DIR_ENVS.joinpath("assessem.yaml")
    resources:
        mem_mb=lambda wildcards, attempt: get_assessem_mem(wildcards.bin_size) * attempt,
        time_hrs=lambda wildcards, attempt: get_assessem_hrs(wildcards.bin_size) * attempt,
    params:
        script=DIR_SCRIPTS.joinpath("asm_label_qc", "assessem-embed.py"),
        bin_size=lambda wildcards: binsize_to_int(wildcards.bin_size),
        fset=lambda wildcards: determine_embed_feature_set(wildcards.feature_set)
    shell:
        "{params.script} --verbose --data-cache {input.cache} "
        "--bin-size {params.bin_size} {params.fset} "
        "--out-binned {output.dataset} --out-trans {output.embed} > {log}"


rule run_all_assessem_jobs:
    # note: a binning resolution of 100 bp
    # can be handled wrt memory (200+ GB but feasible)
    # but requires too much runtime to finish (failed w/ 96 hrs)
    # Hence, drop the 100 bp and compensate by selectively
    # excluding the read depth features from the embedding
    input:
        cache = expand(
            rules.prepare_region_cache.output.hdf,
            sample=SAMPLES
        ),
        embeds = expand(
            rules.compute_embedding.output,
            sample=SAMPLES,
            bin_size=["10k", "1k"],
            feature_set=["full", "no-mq0", "no-rd", "mrgfull"]
        )
