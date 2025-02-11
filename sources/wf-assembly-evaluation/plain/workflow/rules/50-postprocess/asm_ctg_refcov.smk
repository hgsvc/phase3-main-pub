
rule mosdepth_assembly_reference_coverage_window:
    input:
        bam = rules.minimap_assembly_to_reference_align_bam.output.bam,
        bai = rules.minimap_assembly_to_reference_align_bam.output.bai
    output:
        check = DIR_PROC.joinpath(
            "50-postprocess", "asm_ctg_refcov", "mosdepth",
            "{refgenome}", "{sample}.{asm_unit}.mq{mapq}.wd",
            "{sample}.{asm_unit}.{refgenome}.mq{mapq}.ok"
        )
    threads: CPU_LOW
    conda: DIR_ENVS.joinpath("biotools", "mosdepth.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 4096 * attempt,
        time_hrs = lambda wildcards, attempt: max(0, attempt-1)
    params:
        window_size = MOSDEPTH_ASSM_REF_COV_WINDOW_SIZE,
        exclude_flag = 1796,  # same as mosdepth default
        min_mapq = lambda wildcards: int(wildcards.mapq),
        out_prefix = lambda wildcards, output: pathlib.Path(output.check).with_suffix("")
    shell:
        "mosdepth --threads {threads} --by {params.window_size} "
            "--no-per-base --flag {params.exclude_flag} --use-median "
            "--quantize 0:1:2:3:4:5:6: --fast-mode --mapq {params.min_mapq} "
            "{params.out_prefix} {input.bam}"
            " && "
        "touch {output.check}"


rule mosdepth_merge_region_contig_coverage:
    input:
        check = expand(
            rules.mosdepth_assembly_reference_coverage_window.output.check,
            mapq=MOSDEPTH_ASSM_REF_COV_MAPQ_THRESHOLDS,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            allow_missing=True
        )
    output:
        merged_regions = DIR_RES.joinpath(
            "coverage", "contig_ref",
            "{refgenome}", "{sample}.{refgenome}.win-ctg-cov.tsv.gz"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt
    run:
        import pathlib as pl
        import pandas as pd

        merged = []
        for check_file in input.check:
            regions_file = pl.Path(check_file).with_suffix(".regions.bed.gz")
            assert regions_file.is_file(), f"File not found: {regions_file}"
            name_components = regions_file.name.split(".")
            mapq_t = name_components[-4]
            assert mapq_t.startswith("mq")
            mapq_t = mapq_t.upper()
            asm_unit = name_components[-6]
            assert asm_unit.startswith("asm")
            regions = pd.read_csv(
                regions_file, sep="\t", header=None,
                names=["chrom", "start", "end", "ctg_align_cov"],
                index_col=["chrom", "start", "end"]
            )
            regions.columns = pd.MultiIndex.from_tuples(
                [(asm_unit, mapq_t, c) for c in regions.columns],
                names=["asm_unit", "mapq", "statistic"]
            )
            merged.append(regions)
        if len(merged) == 1:
            merged = merged[0]
        else:
            merged = pd.concat(merged, axis=1, ignore_index=False)

        merged.to_csv(output.merged_regions, index=True, header=True, sep="\t")
    # END OF RUN BLOCK


rule run_assembly_reference_coverage:
    """
    TODO
    need to adapt solution with ASSEMBLY_UNITS_NO_CONTAM,
    does not work for special cases such as Verkko's rDNA
    output
    """
    input:
        windowed = expand(
            rules.mosdepth_assembly_reference_coverage_window.output.check,
            refgenome=WILDCARDS_REF_GENOMES,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            mapq=MOSDEPTH_ASSM_REF_COV_MAPQ_THRESHOLDS
        ),
        rdna_win = expand(
            rules.mosdepth_assembly_reference_coverage_window.output.check,
            refgenome=COMPLETE_REF_GENOME,
            sample=SAMPLES,
            asm_unit=ASSEMBLY_UNITS_NO_CONTAM,
            mapq=MOSDEPTH_ASSM_REF_COV_MAPQ_THRESHOLDS
        ),
        merged = expand(
            rules.mosdepth_merge_region_contig_coverage.output.merged_regions,
            refgenome=COMPLETE_REF_GENOME,
            sample=SAMPLES
        )

