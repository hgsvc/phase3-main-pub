
rule mosdepth_assembly_read_coverage_window:
    """Note that this rule assumes a pre-filtered
    BAM as input file. Thus, the exlude flag option
    is explicitly set to 0 (unset).
    """
    input:
        bam = rules.merge_read_to_assembly_subset_alignments.output.bam,
        bai = rules.merge_read_to_assembly_subset_alignments.output.bai
    output:
        check = DIR_PROC.joinpath(
            "50-postprocess", "asm_ctg_readcov", "mosdepth",
            "{sample}.{read_type}.{aln_subset}.mq{mapq}.wd",
            "{sample}.{read_type}.{aln_subset}.mq{mapq}.ok"
        )
    threads: CPU_LOW
    conda: DIR_ENVS.joinpath("biotools", "mosdepth.yaml")
    resources:
        mem_mb = lambda wildcards, attempt: 8192 * attempt,
        time_hrs = lambda wildcards, attempt: max(0, attempt-1)
    params:
        window_size = MOSDEPTH_ASSM_READ_COV_WINDOW_SIZE,
        min_mapq = lambda wildcards: int(wildcards.mapq),
        out_prefix = lambda wildcards, output: pathlib.Path(output.check).with_suffix("")
    shell:
        "mosdepth --threads {threads} --by {params.window_size} "
            "--no-per-base --flag 0 --use-median "
            "--quantize 0:5:8:12:15:20:25:30:50:100: --fast-mode --mapq {params.min_mapq} "
            "{params.out_prefix} {input.bam}"
            " && "
        "touch {output.check}"


rule mosdepth_coverage_stats_summary:
    """ Generate a summary info
    file of the (mean) coverages
    along the assemble sequences
    to be used for normalization
    in later stages

    Makes only sense for primary
    alignments at MAPQ 0
    """
    input:
        check_file = expand(
            rules.mosdepth_assembly_read_coverage_window.output.check,
            aln_subset=["onlyPRI"],
            mapq=["00"],
            allow_missing=True
        ),
        clean_regions = get_clean_assembly_regions
        #clean_regions = rules.define_clean_assembly_regions.output.tag_tig
        # 2024-03-18: changed to realize skipping over contamination scan
    output:
        stats = DIR_PROC.joinpath(
            "50-postprocess", "asm_ctg_readcov", "mosdepth",
            "{sample}.{read_type}.cov-stats.tsv"
        )
    params:
        threshold_percentile = 99.5
    run:
        import pathlib as pl
        import gzip
        import pandas as pd
        import numpy as np
        #import scipy.stats as stats  # Not in default snakemake env
        _this = "50-postprocess::asm_ctg_readcov::summarize_mosdepth_coverage"

        # 2024-03-27 change here to accommodate sequence names that use
        # '#' as a separator - just annoying ...
        if str(input.clean_regions).endswith(".gz"):
            fopen, fmode = gzip.open, "rt"
        else:
            fopen, fmode = open, "r"

        clean_contigs = set()
        with fopen(input.clean_regions, fmode) as listing:
            for line in listing:
                if line.startswith("#") or not line.strip():
                    continue
                contig = line.split()[0]
                clean_contigs.add(contig)

        assert len(input.check_file) == 1  # expand() returns a NamedList
        summary_file = pl.Path(input.check_file[0]).with_suffix(".mosdepth.summary.txt")
        if not summary_file.is_file():
            logerr(f"\nERROR in {_this} - file does not exist {summary_file}\n")
            raise FileNotFoundError(summary_file.name)
        df = pd.read_csv(summary_file, sep="\t", header=0)
        df = df.loc[~df["chrom"].str.endswith("region"), :].copy()
        reported_total_mean = df.loc[df["chrom"] == "total", "mean"].values[0]
        df = df.loc[df["chrom"].isin(clean_contigs), :].copy()  # NB: gets rid of "total" entry

        wt_avg = np.average(df["mean"].values, weights=df["length"].values)
        wt_var = np.average((df["mean"].values - wt_avg)**2, weights=df["length"].values)
        wt_stddev = np.sqrt(wt_var)

        factor = round(params.threshold_percentile / 100, 3)
        assert factor < 1
        threshold_index = int(df["mean"].values.size * factor)
        threshold_score = np.sort(df["mean"].values)[threshold_index]

        assert np.isclose(reported_total_mean, wt_avg, atol=0.5), \
            f"Averages do not match: {reported_total_mean} vs {wt_avg}"

        fmt_pctile = str(params.threshold_percentile).replace(".", "-")
        df = pd.DataFrame(
            [[
                wildcards.sample, wildcards.read_type, reported_total_mean,
                wt_avg, wt_stddev, wt_var, threshold_score
            ]],
            columns=[
                "sample", "read_type", "mosdepth_global_mean_cov",
                "global_mean_cov", "global_mean_cov_stddev",
                "global_mean_cov_var", f"mean_cov_at_pctile_{fmt_pctile}"
            ]
        )
        df.to_csv(output.stats, sep="\t", header=True, index=False)
    # END OF RUN BLOCK


rule transform_mosdepth_window_read_coverage:
    input:
        stats_files = lambda wildcards: expand(
            rules.mosdepth_coverage_stats_summary.output.stats,
            read_type=USE_READ_TYPES_WINDOW_COVERAGE,
            allow_missing=True
        ),
        check_files = lambda wildcards: expand(
            rules.mosdepth_assembly_read_coverage_window.output.check,
            read_type=USE_READ_TYPES_WINDOW_COVERAGE,
            aln_subset=["onlySPL", "onlyPRI"],
            mapq=MOSDEPTH_ASSM_READ_COV_MAPQ_THRESHOLDS,
            allow_missing=True
        ),
        clean_regions = get_clean_assembly_regions
        #clean_regions = rules.define_clean_assembly_regions.output.tag_tig
        # 2024-03-18: changed to realize skipping over contamination scan
    output:
        hdf = DIR_PROC.joinpath(
            "50-postprocess", "asm_ctg_readcov", "mosdepth",
            "{sample}.read-cov.h5"
        )
    resources:
        mem_mb=lambda wildcards, attempt: 4096 * attempt,
        time_hrs=lambda wildcards, attempt: max(0, attempt - 1)
    run:
        import pathlib as pl
        import gzip
        import pandas as pd
        import numpy as np

        # 2024-03-27 change here to accommodate sequence names that use
        # '#' as a separator - just annoying ...
        if str(input.clean_regions).endswith(".gz"):
            fopen, fmode = gzip.open, "rt"
        else:
            fopen, fmode = open, "r"

        clean_contigs = set()
        with fopen(input.clean_regions, fmode) as listing:
            for line in listing:
                if line.startswith("#") or not line.strip():
                    continue
                contig = line.split()[0]
                clean_contigs.add(contig)


        aln_subsets = {"onlyPRI": 1, "onlySPL": 0, "onlySEC": 2}

        global_covs = dict()
        all_stats = []
        for stats_file in sorted(set(input.stats_files)):
            stats = pd.read_csv(stats_file, sep="\t", header=0)
            read_type = pl.Path(stats_file).name.rsplit(".", 3)[-3]
            global_mean_cov = round(stats["global_mean_cov"].values[0], 0)
            global_covs[read_type] = int(global_mean_cov)
            all_stats.append(stats)
        all_stats = pd.concat(all_stats, axis=0, ignore_index=False)

        columns = ["contig", "start", "end", "median_cov"]
        regions = None
        for check_file in sorted(set(input.check_files)):
            region_file = pl.Path(check_file).with_suffix(".regions.bed.gz")
            _, read_type, aln_subset, mapq, _ = pl.Path(check_file).name.rsplit(".", 4)
            aln_type = aln_subsets[aln_subset]
            mapq = int(mapq.strip("mq"))
            new_regions = pd.read_csv(
                region_file, sep="\t", header=None, names=columns,
                index_col=["contig", "start", "end"]
            )
            global_cov = global_covs[read_type]
            new_regions["pct_median_cov"] = (new_regions["median_cov"].values / global_cov * 100).round(1)
            col_idx = pd.MultiIndex.from_tuples(
                [
                    (read_type, aln_type, mapq, "median_cov"),
                    (read_type, aln_type, mapq, "pct_median_cov")
                ],
                names=["read_type", "aln_type", "mapq", "value"]
            )
            new_regions.columns = col_idx
            if regions is None:
                regions = new_regions.copy()
            else:
                regions = regions.join(new_regions, how="outer")

        with pd.HDFStore(output.hdf, "w", complevel=9, complib="blosc") as hdf:
            contig_lut = []
            contig_num = 0
            for contig, windows in regions.groupby("contig"):
                assert isinstance(contig, str)
                if contig not in clean_contigs:
                    continue
                contig_num += 1
                untagged_contig, asm_unit_tag = contig.rsplit(".", 1)
                store_key = f"{asm_unit_tag}/contig{contig_num}"
                hdf.put(store_key, windows, format="fixed")
                contig_lut.append(
                    (contig, untagged_contig, asm_unit_tag, store_key, contig_num)
                )

            contig_lut = pd.DataFrame.from_records(
                contig_lut, columns=[
                    "contig", "untagged_contig", "asm_unit",
                    "store_key", "contig_num"
                ]
            )
            hdf.put("contigs", contig_lut, format="fixed")
            hdf.put("covstats", all_stats, format="fixed")
    # END OF RUN BLOCK


# TODO: need wildcard comb function for
# sample and read type
rule run_all_mosdepth_assembly_read_coverage:
    input:
        check_files = expand(
            rules.mosdepth_assembly_read_coverage_window.output.check,
            sample=SAMPLES,
            read_type=USE_READ_TYPES_WINDOW_COVERAGE,
            aln_subset=["onlyPRI", "onlySPL"],
            mapq=MOSDEPTH_ASSM_READ_COV_MAPQ_THRESHOLDS
        ),
        transform_cov = expand(
            rules.transform_mosdepth_window_read_coverage.output.hdf,
            sample=SAMPLES
        )


rule run_all_mosdepth_coverage_stats:
    input:
        stats = expand(
            rules.mosdepth_coverage_stats_summary.output.stats,
            sample=SAMPLES,
            read_type=USE_READ_TYPES_WINDOW_COVERAGE
        )
